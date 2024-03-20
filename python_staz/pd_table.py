import allel
import pandas as pd
import argparse

# parse parametres
parser = argparse.ArgumentParser(description='Tool to benchmark MGI and ilumina vcfs to GIB true dataset')

parser.add_argument('--true', help='true HG002 variant call file (vcf) downaloaded from GIB') 
parser.add_argument('--mgi', help='MGI vcf')
parser.add_argument('--ilumina', help='ilumina vcf')

args = parser.parse_args()
true=args.true
mgi=args.mgi
ilumina=args.ilumina

# load vcf but missing GT and ALT
print("load vcfs")
callset_true = allel.vcf_to_dataframe(true, fields=['CHROM','POS', 'ID', 'REF', 'ALT', 'GT','QUAL'])
callset_1 = allel.vcf_to_dataframe(mgi, fields=['CHROM','POS', 'ID', 'REF', 'ALT', 'GT','QUAL'])
callset_2 = allel.vcf_to_dataframe(ilumina, fields=['CHROM','POS', 'ID', 'REF', 'ALT', 'GT','QUAL'])
print("load vcfs done")

#df2 = pd.DataFrame(callset)
callset_true.drop(['ALT_2', 'ALT_3'], axis=1, inplace=True)
callset_1.drop(['ALT_2', 'ALT_3'], axis=1, inplace=True)
callset_2.drop(['ALT_2', 'ALT_3'], axis=1, inplace=True)

# rearange GT field
print("load GT")
gt_vcf_true = allel.read_vcf(true, fields = 'GT')
gt_vcf_1 = allel.read_vcf(mgi, fields = 'GT')
gt_vcf_2 = allel.read_vcf(ilumina, fields = 'GT')
print("load GT done")

callset_true.reset_index(drop=True, inplace=True)
callset_1.reset_index(drop=True, inplace=True)
callset_2.reset_index(drop=True, inplace=True)

final_tab_true = pd.concat([callset_true, pd.DataFrame([str(x.flatten().tolist()).replace(', ', '/').replace('[', '').replace(']', '') for x in gt_vcf_true['calldata/GT'][:]])], axis=1)
final_tab_true.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'GT']

final_tab_1 = pd.concat([callset_1, pd.DataFrame([str(x.flatten().tolist()).replace(', ', '/').replace('[', '').replace(']', '') for x in gt_vcf_1['calldata/GT'][:]])], axis=1)
final_tab_1.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'GT']

final_tab_2 = pd.concat([callset_2, pd.DataFrame([str(x.flatten().tolist()).replace(', ', '/').replace('[', '').replace(']', '') for x in gt_vcf_2['calldata/GT'][:]])], axis=1)
final_tab_2.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'GT']


# create and absolute variant identificator columns for all tables (vcfs)
final_tab_true['ident'] = [str(final_tab_true['CHROM'][i] + ":" + str(final_tab_true['POS'][i]) + "_" + final_tab_true['ALT'][i]) for i in range(0, len(final_tab_true['CHROM']))]
final_tab_1['ident'] = [str(final_tab_1['CHROM'][i] + ":" + str(final_tab_1['POS'][i]) + "_" + final_tab_1['ALT'][i]) for i in range(0, len(final_tab_1['CHROM']))]
final_tab_2['ident'] = [str(final_tab_2['CHROM'][i] + ":" + str(final_tab_2['POS'][i]) + "_" + final_tab_2['ALT'][i]) for i in range(0, len(final_tab_2['CHROM']))]

# identification of True/False positive and negative results

## final_table as True dataset
true_pos_1 = sum(final_tab_true['ident'].isin(final_tab_1['ident']))/len(final_tab_true['ident']) * 100
false_pos_1= sum(~final_tab_1['ident'].isin(final_tab_true['ident']))/len(final_tab_1['ident']) * 100
false_neg_1 = sum(~final_tab_true['ident'].isin(final_tab_1['ident']))/len(final_tab_1['ident']) * 100

true_pos_2 = sum(final_tab_true['ident'].isin(final_tab_2['ident']))/len(final_tab_true['ident']) * 100
false_pos_2 = sum(~final_tab_2['ident'].isin(final_tab_true['ident']))/len(final_tab_2['ident']) * 100
false_neg_2 = sum(~final_tab['ident'].isin(final_tab2['ident']))/len(final_tab2['ident']) * 100

print(true_pos_1)
print(false_pos_1)
print(false_neg_1)
