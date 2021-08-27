import re
import sys
import argparse
import edlib

def get_cigar_list(cigarString):
    """
        Returns a list of tuples based on the cigar string. 
        Each tuple has two elements; the first one is showing 
        the operation which could be one of =, X, I or D and the 
        second element is an integer which shows the length of 
        the associated operation.
    """
    cigarOps = re.compile("[0-9]+").split(cigarString)[1:]
    cigarSizes = [int(size) for size in re.compile("=|X|I|D").split(cigarString)[:-1]]
    cigarList = [ (op, size) for op, size in zip(cigarOps, cigarSizes)]
    return cigarList

def is_par(chrom, pos):
    # locations taken from
    # https://uswest.ensembl.org/info/genome/genebuild/human_PARS.html
    X_par1 = (10001, 2781479)
    X_par2 = (155701383, 156030895)
    Y_par1 = (10001, 2781479)
    Y_par2 = (56887903, 57217415)
    if "X" in chrom:
        if pos >= X_par1[0] and pos <= X_par1[1]:
            return True
        if pos >= X_par2[0] and pos <= X_par2[1]:
            return True
    if "Y" in chrom:
        if pos >= Y_par1[0] and pos <= Y_par1[1]:
            return True
        if pos >= Y_par2[0] and pos <= Y_par2[1]:
            return True
    return False

def skip_record(cols):
    if cols[0].startswith("#"):
        return True
    if len(cols) < 9: # skip a rare bug for extremely long variant in dipcall output
        return True
    filter_col = cols[6]
    # Skip 'HET2', 'HET1', 'GAP2', 'GAP1', 'DIPX', 'DIPY' and their combinations
    if filter_col != ".":
        return True
    # skip regions with alignments from only one haplotype
    #(if the vcf is containing only the confident regions
    # it would skip a small number of edge cases)
    if '.' in cols[9]:
        return True

def main():
    parser = argparse.ArgumentParser(description='Estimate Heterozygosity')
    parser.add_argument('--vcf', type=str, help='Path to the vcf file')
    args = parser.parse_args()

    counts_by_region = {'autosome':
                        {'het_snp': 0, 'hom_snp': 0, 'het_snp_ins': 0, 'ins_hap1': 0, 'ins_hap2': 0},
                       'par':
                        {'het_snp': 0, 'hom_snp': 0, 'het_snp_ins': 0, 'ins_hap1': 0, 'ins_hap2': 0},
                       'non_par_X':
                        {'het_snp': 0, 'hom_snp': 0, 'het_snp_ins': 0, 'ins_hap1': 0, 'ins_hap2': 0}}
    with open(args.vcf) as f:
        for line in f:
            cols = line.strip().split()
            if skip_record(cols): # skip some special cases
                continue
            chrom = cols[0]
            pos = int(cols[1])

            # find the corresponding region
            if is_par(chrom, pos):
                region = 'par'
            elif "X" in chrom:
                region = 'non_par_X'
            elif "Y" in chrom:
                continue # skip non-PAR Y
            else:
                region = 'autosome'

            ref = cols[3]
            alt = cols[4].split(",")
            alleles = [ref]
            alleles.extend(alt)
            genotypes = [int(i) for i in cols[9].split(":")[0].split("|")]
            allele_hap1 = "" if alleles[genotypes[0]] == "*" else alleles[genotypes[0]]
            allele_hap2 = "" if alleles[genotypes[1]] == "*" else alleles[genotypes[1]]
            
            if genotypes[0] == genotypes[1]: # homo
                if len(allele_hap1) == 1: # homo snp
                    counts_by_region[region]['hom_snp'] += 1
            elif len(allele_hap1) == 1 and len(allele_hap2) == 1: #snp
                counts_by_region[region]['het_snp'] +=1

            # insertion on both hap1 & hap2 (will do a global alignment between the inserted parts)
            elif len(ref) < len(allele_hap1) and len(ref) < len(allele_hap2):
                result = edlib.align(allele_hap1, allele_hap2, mode="NW", task = "path")
                for op in get_cigar_list(result['cigar']):
                    if op[0] == 'X':
                        counts_by_region[region]['het_snp_ins'] += op[1]
                    elif op[0] == 'I':
                        counts_by_region[region]['ins_hap1'] += op[1]
                    elif op[0] == 'D':
                        counts_by_region[region]['ins_hap2'] += op[1]
            else:
                # deletion on hap1 and insertion (or shorter deletion) on hap2
                if len(allele_hap1) < len(allele_hap2):
                    counts_by_region[region]['ins_hap2'] += len(allele_hap2) - len(allele_hap1)
                # deletion on hap2 and insertion (or shorter deletion) on hap1
                else:
                    counts_by_region[region]['ins_hap1'] += len(allele_hap1) - len(allele_hap2)
    
    # print results
    print("#{}\t{}\t{}\t{}\t{}\t{}".format('region', 'het_snp', 'hom_non_ref_snp', 'het_snp_within_insertion', 'insertion_hap1', 'ins_hap2'))
    for region in counts_by_region:
        c = counts_by_region[region]
        print("{}\t{}\t{}\t{}\t{}\t{}".format(region, c['het_snp'], c['hom_snp'], c['het_snp_ins'], c['ins_hap1'], c['ins_hap2']))
    het_ratio = (counts_by_region['autosome']['het_snp'] +  counts_by_region['par']['het_snp']) / (counts_by_region['autosome']['hom_snp'] + counts_by_region['par']['hom_snp'] )
    print("\nHet Ratio (het_snp / hom_non_ref_snp) in autosome and par: {:.2f}".format(het_ratio))


if __name__ == "__main__":
    main()
