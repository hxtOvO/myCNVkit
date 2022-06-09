"""DataFrame-level subdivide operation.

Split each region into similar-sized sub-regions.

The functions here operate on pandas DataFrame and Series instances, not
GenomicArray types.

"""
import logging
import random

import pandas as pd

from .merge import merge


def subdivide(table, avg_size, min_size=0, verbose=False):
    return pd.DataFrame.from_records(
        _split_targets(table, avg_size, min_size, verbose),
        columns=table.columns)


def _split_targets(regions, avg_size, min_size, verbose):
    """Split large regions into smaller, consecutive regions.

    Output bin metadata and additional columns match the input dataframe.

    Parameters
    ----------
    avg_size : int
        Split regions into equal-sized subregions of about this size.
        Specifically, subregions are no larger than 150% of this size, no
        smaller than 75% this size, and the average will approach this size when
        subdividing a large region.
    min_size : int
        Drop any regions smaller than this size.
    verbose : bool
        Print a log message when subdividing a region.

    """
    for row in merge(regions).itertuples(index=False):  # Merge overlapping rows in a DataFrame 迭代器
        span = row.end - row.start
        if span >= min_size:  # 如果区间bin长度大于最小，那么需要将其划分为多个bins
            nbins = int(round(span / avg_size)) or 1
            if nbins == 1:
                yield row
            else:
                # Divide the region into equal-sized bins 这就是固定区间的划分策略，需要将其修改为贴近变异区间的位置
                bin_size = span / nbins
                bin_start = row.start
                if verbose:  # 如果需要同时输出log信息
                    label = (row.gene if 'gene' in regions else
                             "%s:%d-%d" % (row.chromosome, row.start, row.end))
                    logging.info("Splitting: {:30} {:7} / {} = {:.2f}"
                                 .format(label, span, nbins, bin_size))
                for i in range(1, nbins):
                    bin_end = row.start + int(i * bin_size)
                    yield row._replace(start=bin_start, end=bin_end)
                    bin_start = bin_end
                yield row._replace(start=bin_start)


                #  Divide the region into random-sized bins
                #  ignore its nbins, random split
                # binlist = _random_split(row.start, row.end, avg_size)
                # bin_start = row.start
                # if verbose:  # 如果需要同时输出log信息
                #     label = (row.gene if 'gene' in regions else
                #              "%s:%d-%d" % (row.chromosome, row.start, row.end))
                #     logging.info("Splitting: {:30} {:7} / {} "
                #                  .format(label, span, nbins))
                # for i in range(1, binlist.__len__()):
                #     # print(binlist[i - 1])
                #     bin_end = bin_start + int(binlist[i - 1])
                #     yield row._replace(start=bin_start, end=bin_end)
                #     bin_start = bin_end
                # yield row._replace(start=bin_start)


def _random_split(begin, end, avg_size):
    """Split large regions into smaller, consecutive regions.

    Output: An array with binsize of every bin

    Parameters
    ----------
    begin : int
        row.begin
    end : int
        row.end
    """
    span = end - begin
    list = []
    for i in range(1, int(span)):
        len = int(random.randint(int(avg_size*0.95),int(avg_size*1.05)))
        span = span - len
        if span >= 0:
            list.append(len)
        else:
            list.append(span+len)
            break
    # logging.info(list)
    return list
