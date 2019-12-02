"""
Clean and Extract chromosome 20 variant calls only from Illumina, PacBio and Truth VCF.
Store them into CSV for further process
_experiment_= "1"
__author__ = “Jiajia Xu”
__copyright__ = “Copyright 2019, Individual Computing Project”
__email__ = “jiajia.xu@anu.edu.au”
__status__ = “Production”
"""

# This function cleans and extracts chromosome 20 from Illumina and pacbio VCF
def HG_cut_chr20(filename):
    file = open(filename, "r")
    variants = []
    variant=[]
    count = 0
    pos_list=[]
    for line in file:
        if line.startswith('#'):
            continue

        rows = line.split()
        chrom = rows[0].strip()
        pos = int(rows[1])
        ID=rows[2].strip()
        ref = rows[3].strip()
        alt = rows[4].strip()
        qual=rows[5].strip()
        filter=rows[6].strip()
        info=rows[7].strip()
        format=rows[8].strip().split(':')
        HG002=rows[9].strip().split(':')

        if chrom!='20':
            continue

        assert len(format)==6
        genotype=HG002[format.index('GT')]
        quality=HG002[format.index('GQ')]
        depth=HG002[format.index('DP')]
        depth_each_allele=HG002[format.index('AD')]
        variant_allele_fractions=HG002[format.index('VAF')]
        phred_scaled_likelihoods=HG002[format.index('PL')]

        variant = [chrom, pos, ref, alt, info, qual,filter,info,genotype,quality,depth,depth_each_allele,
                   variant_allele_fractions,phred_scaled_likelihoods]

        # print(variant)
        variants.append(variant)
        pos_list.append(pos)
        count = count + 1
    file.close()
    return pos_list, count, variants


# This function cleans and extracts chromosome 20 from Truth VCF
def truth_cut_chr20(filename):
    file = open(filename, "r")
    variants = []
    variant=[]
    count = 0
    pos_list=[]
    for line in file:
        if line.startswith('#'):
            continue

        rows = line.split()
        chrom = rows[0].strip()
        pos = int(rows[1])
        ID=rows[2].strip()
        ref = rows[3].strip()
        alt = rows[4].strip()
        qual=rows[5].strip()
        filter=rows[6].strip()
        info=rows[7].strip()
        format=rows[8].strip().split(':')
        HG002=rows[9].strip().split(':')

        if chrom!='20':
            continue

        assert len(format)==8
        genotype=HG002[format.index('GT')]
        depth = HG002[format.index('DP')]
        depth_all_allele = HG002[format.index('ADALL')]
        depth_each_allele = HG002[format.index('AD')]
        quality=HG002[format.index('GQ')]
        origin_genotype = HG002[format.index('IGT')]
        phase_IGT = HG002[format.index('IPS')]
        phase_GT = HG002[format.index('PS')]

        variant = [chrom, pos, ref, alt, info, qual, filter, info, genotype, depth,depth_all_allele,
                   depth_each_allele, quality, origin_genotype, phase_IGT,phase_GT]

        # print(variant)
        variants.append(variant)
        pos_list.append(pos)
        count = count + 1

    file.close()
    return pos_list, count, variants


## Use the above functions to get chromosome 20 variant calls from VCFs (Illumina, Pacbio and Truth)
p1,count1,wgs_chr20_out=HG_cut_chr20('../../raw_Data/wgs.HG002.output.vcf')
p2,count2,pacbio_chr20_out=HG_cut_chr20('../../raw_Data/pacbio.HG002.output.vcf')
p3,count3,truth_chr20=truth_cut_chr20('../../raw_Data/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf')



#Store chromosome 20 variant calls as CSV
import csv

with open("wgs.chr20.output.csv","w") as file:
    writer=csv.writer(file)
    writer.writerow(['chrom','pos','ref','alt','info','qual','filter','info','genotype','quality','depth',
                     'depth_each_allele','variant_allele_fractions','phred_scaled_likelihoods'])
    writer.writerows(wgs_chr20_out)

with open("pacbio.chr20.output.csv","w") as file:
    writer=csv.writer(file)
    writer.writerow(['chrom','pos','ref','alt','info','qual','filter','info','genotype','quality','depth',
                     'depth_each_allele','variant_allele_fractions','phred_scaled_likelihoods'])
    writer.writerows(pacbio_chr20_out)

with open("truth.chr20.csv","w") as file:
    writer=csv.writer(file)
    writer.writerow(['chrom','pos','ref','alt','info','qual', 'filter', 'info', 'genotype', 'depth','depth_all_allele',
                   'depth_each_allele', 'quality', 'origin_genotype', 'phase_IGT','phase_GT'])
    writer.writerows(truth_chr20)
