import allel
import pandas as pd

# Načtení dat z VCF souboru
callset = allel.read_vcf("HG002_HiSeq_sorted_mapped_haplCall.g.vcf.gz")

# Extrahování informací o genotypu vzorků
genotypes = allel.GenotypeArray(callset['calldata/GT'])

# Vytvoření tabulky pro genotypy
genotype_table = pd.DataFrame(genotypes, columns=callset['samples'])
genotype_table.index.name = 'Variant'

# Připojení tabulky genotypů k původní tabulce
table_with_genotypes = pd.concat([table, genotype_table], axis=1)

# Zobrazení prvních pěti řádků
print(table_with_genotypes.head())
