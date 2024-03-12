import allel
import pandas as pd

# Načtení dat z VCF souboru
callset = allel.read_vcf("HG002_HiSeq_sorted_mapped_haplCall.g.vcf.gz")

# Extrahování informací o genotypu
# GenotypeArray - dat. struktura z knihovny scikit-allel, umožňuje ukládat a manipulovat z genotypy z .vcf souboru 
genotypes = allel.GenotypeArray(callset['calldata/GT'])

# Vytvoření tabulky pro genotypy
gt_table = pd.DataFrame(genotypes, columns=callset['samples'])
gt_table.index.name = 'Variant'

# Připojení tabulky genotypů k původní tabulce
table_with_gt = pd.concat([table, gt_table], axis=1)

print(table_with_gt.head())
