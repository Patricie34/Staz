import allel

# Because working with genotype calls is a very common task, scikit-allel has a GenotypeArray class 
# which adds some convenient functionality to an array of genotype calls. 
# To use this class, pass the raw NumPy array into the GenotypeArray class constructor, e.g.:
callset = allel.read_vcf('example.vcf')
gt = allel.GenotypeArray(callset['calldata/GT'])
print(gt)
gt.is_het()

# the is_het() method locates all heterozygous genotype calls
print(gt.is_het())  

# the count_het() method will count heterozygous calls, summing over variants (axis=0) or samples (axis=1) if requested. 
# E.g., to count the number of het calls per variant:
print(gt.count_het(axis=1)) 

# to extract absolutely everything from a VCF file, then you can provide a special value '*' as the fields parameter:
callset = allel.read_vcf('example.vcf', fields='*')
print(sorted(callset.keys()))
