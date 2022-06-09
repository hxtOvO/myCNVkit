"""Supporting functions for the 'antitarget' command."""
import logging
import math
import os.path
import time
from concurrent import futures
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
from io import StringIO
from skgenome import tabio

from . import core, samutil
from .cnary import CopyNumArray as CNA
from .parallel import rm, to_chunks
from .params import NULL_LOG2_COVERAGE

plot = []

def do_coverage(bed_fname, bam_fname, by_count=False, min_mapq=0, processes=1):
    """Calculate coverage in the given regions from BAM read depths."""
    if not samutil.ensure_bam_sorted(bam_fname):
        raise RuntimeError("BAM file %s must be sorted by coordinates"
                           % bam_fname)
    samutil.ensure_bam_index(bam_fname)
    # ENH: count importers.TOO_MANY_NO_COVERAGE & warn
    cnarr = interval_coverages(bed_fname, bam_fname, by_count, min_mapq,
                               processes)

    return cnarr


def interval_coverages(bed_fname, bam_fname, by_count, min_mapq, processes):
    """Calculate log2 coverages in the BAM file at each interval."""
    meta = {'sample_id': core.fbase(bam_fname)}
    start_time = time.time()

    # Skip processing if the BED file is empty
    with open(bed_fname) as bed_handle:
        for line in bed_handle:
            if line.strip():
                break
        else:
            logging.info("Skip processing %s with empty regions file %s",
                         os.path.basename(bam_fname), bed_fname)
            return CNA.from_rows([], meta_dict=meta)

    # Calculate average read depth in each bin
    if by_count:
        results = interval_coverages_count(bed_fname, bam_fname, min_mapq,
                                           processes)
        read_counts, cna_rows = zip(*results)
        read_counts = pd.Series(read_counts)
        cnarr = CNA.from_rows(list(cna_rows),
                              columns=CNA._required_columns + ('depth',),
                              meta_dict=meta)
    elif "normal" in bam_fname:
        #  如果是一个normal样本，则不对其进行启发式区间处理
        table = interval_coverages_pileup(bed_fname, bam_fname, min_mapq,
                                          processes)
        read_len = samutil.get_read_length(bam_fname)
        read_counts = table['basecount'] / read_len
        table = table.drop('basecount', axis=1)

        #  print table
        f = open("table " + bam_fname + ".txt", "w")
        f.write(str(table))
        f.close()

        cnarr = CNA(table, meta)
        tabio.write(cnarr, 'normal_debug.cnr')
    else:
        #  对tumor样本进行启发式区间划分合并策略
        table = interval_coverages_pileup(bed_fname, bam_fname, min_mapq,
                                          processes)
        read_len = samutil.get_read_length(bam_fname)
        read_counts = table['basecount'] / read_len
        table = table.drop('basecount', axis=1)
        #  print table
        f = open("table " + bam_fname + ".txt", "w")
        f.write(str(table))
        f.close()
        plot.append(len(table))

        #  在这里进行区间合并
        # 根据depth进行启发式的合并
        # 考虑一个典型的启发式算法“模拟退火算法”。
        # 设当前有的区间是A，A要合并了相邻区间B，
        # 情况1：合并后A+B的平均深度depth增加，说明A+B中更可能包含cnv，此时，算法执行合并操作；
        # 否则，情况2：合并后A+B的平均depth降低，说明A+B中也可能包含cnv（此时合并B后的下降是由于噪声和随机性），也可能不包含
        # 此时，根据下降的比例算出一个合并概率，平均深度下降的越多，合并的概率越低。

        #  先写一个简单的两路合并
        # tmp = merge_two_bins(table)

        #  这里生成一个合并相邻区间方案
        tmp_method = p_divide_method(table)
        #  根据合并方案执行合并
        tmp = p_divide_action(table, tmp_method)
        length = len(tmp_method)

        #  得到了初始解后开始迭代
        while True:
            cmp = tmp_method  # cmp是上一步的绝对划分方案
            tmp_method = p_divide_method(tmp)  # tmp_method是正进行这步的相对于上一个状态的划分
            #  关键是在这里把相对method转化为绝对method
            abs_method = to_abs(cmp, tmp_method)
            #  如果此时的method已经平衡或者个数过少，则停止退火
            if len(tmp_method) == length or len(tmp_method) <= 50:
                tmp_method = cmp  # 退出循环时method得记录上一步的绝对method才行
                break

            tmp = p_divide_action(tmp, tmp_method)
            tmp_method = abs_method  # 使用相对method产生新的table后，把这一步的绝对method传递给下一个循环
            length = len(tmp_method)

        #  开始绘制马尔科夫链
        x = [i for i in range(0, len(plot))]
        # y_2 = [compare for i in range(0, len(plot))]
        plt.plot(x, plot, marker='o',color='r', label='Modified')
        # plt.plot(x, y_2, marker='*',color='b', label='标准CNVkit')
        plt.legend()
        plt.xlabel('Number of iterations')
        plt.ylabel('bins count')
        plt.savefig('markf.png')

        f = open("method.txt", "w")
        f.write(str(tmp_method))
        f.close()

        f = open("tmp " + bam_fname + ".txt", "w")
        f.write(str(tmp))
        f.close()

        cnarr = CNA(tmp, meta)
        tabio.write(cnarr, 'tumor_debug.cnr')

    # Log some stats
    tot_time = time.time() - start_time
    tot_reads = read_counts.sum()
    logging.info("Time: %.3f seconds (%d reads/sec, %s bins/sec)",
                 tot_time,
                 int(round(tot_reads / tot_time, 0)),
                 int(round(len(read_counts) / tot_time, 0)))
    logging.info("Summary: #bins=%d, #reads=%d, "
                 "mean=%.4f, min=%s, max=%s ",
                 len(read_counts),
                 tot_reads,
                 (tot_reads / len(read_counts)),
                 read_counts.min(),
                 read_counts.max())
    tot_mapped_reads = samutil.bam_total_reads(bam_fname)
    if tot_mapped_reads:
        logging.info("Percent reads in regions: %.3f (of %d mapped)",
                     100. * tot_reads / tot_mapped_reads,
                     tot_mapped_reads)
    else:
        logging.info("(Couldn't calculate total number of mapped reads)")
    # print(str(type(cnarr)))
    return cnarr


def to_abs(last_method, method):
    """TODO:将相对method，参考abs_method转化为绝对method"""
    abs_method = []
    logging.info("method max = %s , and last_method len is %s", method[len(method)-1][1], len(last_method))
    global plot
    plot.append(len(last_method))
    for i in range(0, len(method)):
        start = method[i][0]
        end = method[i][1]
        abs_method.append([last_method[start][0], last_method[end][1]])
    return abs_method


def p_divide_method(table):
    """"TODO: 生成一个划分方案，一个二维数组：[[0,2],[3,4],[5,5]]
            表示将0,1,2行划分为一个区间，3，4划分为一个，5行自己作为一个区间
    """
    logging.info("processing p divide method ..")
    method = []

    def divide_merge(bin_table, start):
        """"判断从start行开始合并到哪一行，返回合并的结束行"""
        end = start + 1
        if start >= len(bin_table) or end >= len(bin_table):
            return start
        is_accept = True
        count = 1
        while judge_two_bins(bin_table, start, end):
            #  如果合并成功
            end += 1
            if count > 0:
                start += 1
                count = count-1
            if start >= len(bin_table) or end >= len(bin_table):
                break
        # logging.info("this time end is %d", end-1)
        return end - 1

    def judge_two_bins(bin_table, start, end):
        """"根据table里的depth信息判断start行与end行的合并情况与合并概率，返回是否合并bool类型"""
        depth = bin_table.loc[start]['depth']
        next_depth = bin_table.loc[end]['depth']
        # distance = 0.0
        if depth == 0 or next_depth == 0:
            #  如果出现 depth = 0，这时不能直接做除法
            distance = abs(depth - next_depth)
            if distance <= 1.5:
                return True
            else:
                # np.random.seed(0)
                #  p表示接受合并的概率，这个概率应该与depth变化程度，区间跨度相关
                #  depth变化越小，区间跨度越小，接受概率p越大；反之，depth变化大，区间跨度越大，越不容易接受
                p = math.exp(-distance * (end - start))
                p_box = np.array([p, 1 - p])
                index = np.random.choice([True, False], p=p_box.ravel())
                return index
        else:
            #  计算合并后的depth和log2
            span1 = bin_table.loc[start]['end'] - bin_table.loc[start]['start']
            span2 = bin_table.loc[end]['end'] - bin_table.loc[end]['start']
            new_depth = (bin_table.loc[start]['depth'] * span1 + bin_table.loc[end]['depth'] * span2) / (span1 + span2)
            new_log2 = np.log2(new_depth)
            #  如果合并后的new_depth与原来的depth接近(在0.95到1.05倍之间)，则接受合并区间
            #  否则，以p=abs(new_depth-depth)/depth*10 接受合并区间，1-p拒绝合并区间
            distance = abs(new_depth / depth)

            if 0.99 <= distance <= 1.01:
                return True
            if distance < 0.75 or distance > 1.25:
                return False
            else:
                # np.random.seed(0)
                #  p表示接受合并的概率，这个概率应该与depth变化程度，区间跨度相关
                #  depth变化越小，区间跨度越小，接受概率p越大；反之，depth变化大，区间跨度越大，越不容易接受
                #  e^(1/2)*4*(x-0.75)*e^(-1/2*16*(x-0.75)^2)

                p = 4*math.exp(0.5)*(distance-0.75) * math.exp((-(distance-0.75) * (distance-0.75)) * 8)
                if p < 0:
                    return False
                else:
                    p_box = np.array([p, 1 - p])
                    index = np.random.choice([True, False], p=p_box)
                    # index = False
                    # logging.info("the p here is %.3f and choose %s", p, str(index))
                    return index

    bin_start = 0
    while True:
        bin_end = divide_merge(table, bin_start)
        method.append([bin_start, bin_end])
        bin_start = bin_end + 1
        # logging.info("merge %d", bin_start)
        if bin_start >= len(table):
            logging.info("method finished.")
            break
    return method


def p_divide_action(table, method):
    """"TODO: 根据合并方案method合并区间，返回合并后的table"""
    logging.info("handling p divide method ...")
    tmp = pd.DataFrame(columns=table.columns)
    for i in range(0, len(method)):
        #  遍历method中的所有区间合并方案，把每个区间执行合并到tmp中
        start = method[i][0]
        end = method[i][1]
        depth_sum = 0
        for i in range(start, end + 1):
            depth_sum += table.loc[i]['depth'] * (table.loc[i]['end'] - table.loc[i]['start'])
            # span_sum += table.loc[i]['end'] - table.loc[i]['start']
        new_depth = depth_sum / (table.loc[end]['end'] - table.loc[start]['start'])
        if new_depth == 0:
            cmp = 2**table.loc[start]['log2'] + 2**table.loc[end]['log2']
            new_log2 = np.log2(cmp)
        else:
            new_log2 = np.log2(new_depth)
        list = pd.DataFrame([[table.loc[start]['chromosome'], int(table.loc[start]['start']),
                              int(table.loc[end]['end']), '-', new_depth, new_log2]],
                            columns=tmp.columns)
        tmp = tmp.append(list, ignore_index=True)
    # logging.info(tmp.dtypes)
    #  规范dataframe的数据结构
    tmp['chromosome'] = tmp['chromosome'].astype(np.str)
    tmp['start'] = tmp['start'].astype(np.int)
    tmp['end'] = tmp['end'].astype(np.int)
    tmp['depth'] = tmp['depth'].astype(np.float)
    tmp['log2'] = tmp['log2'].astype(np.float)
    # logging.info(tmp.dtypes)
    return tmp


def merge_two_bins(table):
    logging.info("merging two bins ..")
    tmp = pd.DataFrame(columns=table.columns)
    for i in range(0, len(table), 2):
        span1 = table.loc[i]['end'] - table.loc[i]['start']
        span2 = table.loc[i + 1]['end'] - table.loc[i + 1]['start']
        depth = (table.loc[i]['depth'] * span1 + table.loc[i + 1]['depth'] * span2) / (span1 + span2)
        log2 = np.log2(depth)
        list = pd.DataFrame([['chr19', int(table.loc[i]['start']), int(table.loc[i + 1]['end']), '-', depth, log2]],
                            columns=tmp.columns)
        # print(str(list))
        tmp = tmp.append(list, ignore_index=True)
    # logging.info(tmp.dtypes)
    tmp['start'] = tmp['start'].astype(np.int)
    tmp['end'] = tmp['end'].astype(np.int)
    tmp['depth'] = tmp['depth'].astype(np.float)
    tmp['log2'] = tmp['log2'].astype(np.float)
    # logging.info(tmp.dtypes)
    return tmp


def interval_coverages_count(bed_fname, bam_fname, min_mapq, procs=1):
    """Calculate log2 coverages in the BAM file at each interval."""
    regions = tabio.read_auto(bed_fname)
    if procs == 1:
        bamfile = pysam.Samfile(bam_fname, 'rb')
        for chrom, subregions in regions.by_chromosome():
            logging.info("Processing chromosome %s of %s",
                         chrom, os.path.basename(bam_fname))
            for count, row in _rdc_chunk(bamfile, subregions, min_mapq):
                yield [count, row]
    else:
        with futures.ProcessPoolExecutor(procs) as pool:
            args_iter = ((bam_fname, subr, min_mapq)
                         for _c, subr in regions.by_chromosome())
            for chunk in pool.map(_rdc, args_iter):
                for count, row in chunk:
                    yield [count, row]


def _rdc(args):
    """Wrapper for parallel."""
    return list(_rdc_chunk(*args))


def _rdc_chunk(bamfile, regions, min_mapq):
    if isinstance(bamfile, str):
        bamfile = pysam.Samfile(bamfile, 'rb')
    for chrom, start, end, gene in regions.coords(["gene"]):
        yield region_depth_count(bamfile, chrom, start, end, gene, min_mapq)


def region_depth_count(bamfile, chrom, start, end, gene, min_mapq):
    """Calculate depth of a region via pysam count.

    i.e. counting the number of read starts in a region, then scaling for read
    length and region width to estimate depth.

    Coordinates are 0-based, per pysam.
    """

    def filter_read(read):
        """True if the given read should be counted towards coverage."""
        return not (read.is_duplicate
                    or read.is_secondary
                    or read.is_unmapped
                    or read.is_qcfail
                    or read.mapq < min_mapq)

    count = 0
    bases = 0
    for read in bamfile.fetch(reference=chrom, start=start, end=end):
        if filter_read(read):
            count += 1
            # Only count the bases aligned to the region
            rlen = read.query_alignment_length
            if read.pos < start:
                rlen -= start - read.pos
            if read.pos + read.query_alignment_length > end:
                rlen -= read.pos + read.query_alignment_length - end
            bases += rlen
    depth = bases / (end - start) if end > start else 0
    row = (chrom, start, end, gene,
           math.log(depth, 2) if depth else NULL_LOG2_COVERAGE,
           depth)
    return count, row


def interval_coverages_pileup(bed_fname, bam_fname, min_mapq, procs=1):
    """Calculate log2 coverages in the BAM file at each interval."""
    logging.info("Processing reads in %s", os.path.basename(bam_fname))
    if procs == 1:
        table = bedcov(bed_fname, bam_fname, min_mapq)
    else:
        chunks = []
        with futures.ProcessPoolExecutor(procs) as pool:
            args_iter = ((bed_chunk, bam_fname, min_mapq)
                         for bed_chunk in to_chunks(bed_fname))
            for bed_chunk_fname, table in pool.map(_bedcov, args_iter):
                chunks.append(table)
                rm(bed_chunk_fname)
        table = pd.concat(chunks, ignore_index=True)
    # Fill in CNA required columns
    if 'gene' in table:
        table['gene'] = table['gene'].fillna('-')
    else:
        table['gene'] = '-'
    # User-supplied bins might be zero-width or reversed -- skip those
    spans = table.end - table.start
    ok_idx = (spans > 0)
    table = table.assign(depth=0, log2=NULL_LOG2_COVERAGE)
    table.loc[ok_idx, 'depth'] = (table.loc[ok_idx, 'basecount']
                                  / spans[ok_idx])
    ok_idx = (table['depth'] > 0)
    table.loc[ok_idx, 'log2'] = np.log2(table.loc[ok_idx, 'depth'])
    return table


def _bedcov(args):
    """Wrapper for parallel."""
    bed_fname = args[0]
    table = bedcov(*args)
    return bed_fname, table


def bedcov(bed_fname, bam_fname, min_mapq):
    """Calculate depth of all regions in a BED file via samtools (pysam) bedcov.

    i.e. mean pileup depth across each region.
    """
    # Count bases in each region; exclude low-MAPQ reads
    cmd = [bed_fname, bam_fname]
    if min_mapq and min_mapq > 0:
        cmd.extend(['-Q', bytes(min_mapq)])
    try:
        raw = pysam.bedcov(*cmd, split_lines=False)
    except pysam.SamtoolsError as exc:
        raise ValueError("Failed processing %r coverages in %r regions. "
                         "PySAM error: %s" % (bam_fname, bed_fname, exc))
    if not raw:
        raise ValueError("BED file %r chromosome names don't match any in "
                         "BAM file %r" % (bed_fname, bam_fname))
    columns = detect_bedcov_columns(raw)
    table = pd.read_csv(StringIO(raw), sep='\t', names=columns, usecols=columns)
    return table


def detect_bedcov_columns(text):
    """Determine which 'bedcov' output columns to keep.

    Format is the input BED plus a final appended column with the count of
    basepairs mapped within each row's region. The input BED might have 3
    columns (regions without names), 4 (named regions), or more (arbitrary
    columns after 'gene').
    """
    firstline = text[:text.index('\n')]
    tabcount = firstline.count('\t')
    if tabcount < 3:
        raise RuntimeError("Bad line from bedcov:\n%r" % firstline)
    if tabcount == 3:
        return ['chromosome', 'start', 'end', 'basecount']
    if tabcount == 4:
        return ['chromosome', 'start', 'end', 'gene', 'basecount']
    # Input BED has arbitrary columns after 'gene' -- ignore them
    fillers = ["_%d" % i for i in range(1, tabcount - 3)]
    return ['chromosome', 'start', 'end', 'gene'] + fillers + ['basecount']
