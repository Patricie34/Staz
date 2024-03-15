import allel
import pandas as pd

# specify directory where vcf occur -> specify only if using linux container
dir="/home/user/project/results/test/vcf_test/"

# load vcf but missing GT and ALT
callset = allel.vcf_to_dataframe(dir + "TEST_HaplotypeCaller_HG002_HiSeq_2MergedRuns.vcf", fields=['CHROM','POS', 'ID', 'REF', 'ALT', 'GT','QUAL'])


#df2 = pd.DataFrame(callset)
callset.drop(['ALT_2', 'ALT_3'], axis=1, inplace=True)

# rearange GT field
gt_vcf = allel.read_vcf(dir + "TEST_HaplotypeCaller_HG002_HiSeq_2MergedRuns.vcf", fields = 'GT')

final_tab = pd.concat([callset, pd.DataFrame([str(x.flatten().tolist()).replace(', ', '/').replace('[', '').replace(']', '') for x in gt_vcf['calldata/GT'][:]])], axis=1)
final_tab.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'GT']

# add column 'VAR_ID'
>>> test['CHROM'][0] + ":" + str(test['POS'][0]) + "_" + test['ALT_1'][0]
# 'chr1:10492_T'

# print(callset.keys())
# dict_keys = (['samples', 'calldata/GT', 'variants/ALT', 'variants/CHROM', 'variants/FILTER_PASS', 'variants/ID', 'variants/POS', 'variants/QUAL', 
# 'variants/REF'])

# table = allel.vcf_to_dataframe("HG002_HiSeq_sorted_mapped_haplCall.g.vcf.gz")

# Tabulka variant
#table = allel.vcf_to_dataframe("HG002_HiSeq_sorted_mapped_haplCall.g.vcf.gz", fields=['CHROM','POS', 'ID', 'REF', 'ALT_1', 'ALT_2', 'ALT_3', 'QUAL','FILTER_PASS'])
#table.set_index('CHROM', inplace=True)

print(table.head())
#print(table.columns)
