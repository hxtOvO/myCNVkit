from skgenome import tabio
import logging


def merge(table):
    length = len(table)
    table.append([0, 0])
    tmp = []
    i = 0
    while i < length:
        j = i+1
        flag = True
        while j <= length:
            if table[j-1][1] == table[j][0]:
                j = j+1
                continue
            else:
                tmp.append([table[i][0], table[j-1][1]])
                i = j
                flag = False
                break
        if flag:
            i = i+1
        else:
            continue
    return tmp


print("processing compare batch on results and standard results...")
cnarr = tabio.read("results/tumor_sort.cnr")
standard_cnarr = tabio.read("stdard_results/tumor_sort.cnr")

repeat = []
delete = []
standard_repeat = []
standard_delete = []
# print("len of results: %d", len(cnarr))
# print("len of standard: %d", len(standard_cnarr))
for i in range(0, len(cnarr)):
    if cnarr[i]['log2'] >= 1:
        # 重复
        repeat.append([cnarr[i]['start'], cnarr[i]['end']])
    elif cnarr[i]['log2'] <= -1:
        # 缺失
        delete.append([cnarr[i]['start'], cnarr[i]['end']])
    else:  # 正常
        continue
for i in range(0, len(standard_cnarr)):
    if standard_cnarr[i]['log2'] >= 1:
        # 重复
        standard_repeat.append([standard_cnarr[i]['start'], standard_cnarr[i]['end']])
    elif standard_cnarr[i]['log2'] <= -1:
        # 缺失
        standard_delete.append([standard_cnarr[i]['start'], standard_cnarr[i]['end']])
    else:  # 正常
        continue

repeat = merge(repeat)
delete = merge(delete)
standard_repeat = merge(standard_repeat)
standard_delete = merge(standard_delete)

for i in range(0, len(repeat)):
    print("detect CNV duplication in bins:", repeat[i][0], repeat[i][1])
for i in range(0, len(delete)):
    print("detect CNV deletion in bins:", delete[i][0], delete[i][1])
for i in range(0, len(standard_repeat)):
    print("standardCNVkit detect CNV duplication in bins:", standard_repeat[i][0], standard_repeat[i][1])
for i in range(0, len(standard_delete)):
    print("standardCNVkit detect CNV deletion in bins:", standard_delete[i][0], standard_delete[i][1])
