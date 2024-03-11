import allel
import pandas as pd

callset = allel.read_vcf("HG002_HiSeq_sorted_mapped_haplCall.g.vcf.gz")
# print(callset.keys())
# dict_keys = (['samples', 'calldata/GT', 'variants/ALT', 'variants/CHROM', 'variants/FILTER_PASS', 'variants/ID', 'variants/POS', 'variants/QUAL', 
# 'variants/REF'])

# table = allel.vcf_to_dataframe("HG002_HiSeq_sorted_mapped_haplCall.g.vcf.gz")

table = allel.vcf_to_dataframe("HG002_HiSeq_sorted_mapped_haplCall.g.vcf.gz", fields=['CHROM','POS', 'ID', 'REF', 'ALT_1', 'ALT_2', 'ALT_3', 'QUAL','FILTER_PASS'])
table.set_index('CHROM', inplace=True)

print(table.head())
#print(table.columns)
