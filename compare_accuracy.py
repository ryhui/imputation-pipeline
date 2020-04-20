import numpy as np
import pandas as pd
import os
import subprocess

chrom_list = np.arange(1, 23)
mafbins = [0.001, 0.01, 0.05, 0.1, 0.3]
truth_vcf = "full_call/NE1.chr{}.hc.default.gq30_dp10.vcf" # VCF containing genotypes called from the original sequences 
maf_file = "1kg_af/chr{}.af.txt" # Allele frequencies in 1KG, columns: CHROM, POS, AF

cov = [0.05, 0.1, 0.5, 0.75, 1, 1.5, 2]
test_vcf_list = ["NE1.chr{{}}.atlas.{cov}X_MaximumLikelihood.imputed.snp.vcf".format(cov = c) for c in cov] # List of VCFs containing imputed genotypes

sample_name = test_vcf_list[0].split('.')[0] # Assuming that all the VCFs to be tested has the same prefix before '.', which is taken as the sample name by snpSift

def get_stats(res, maf):
    res['maf'] = pd.Series(np.array([maf[s] if s in maf.keys() else np.nan for s in res.index.values]), index = res.index)
    res = res[res['maf'].notnull()]
    res['maf_group'] = np.digitize(res['maf'], mafbins, right = True)
    res = res[res['maf_group'] > 0]
    counts = res.groupby('maf_group').sum()
    counts.sort_index(inplace = True)
    het_correct = counts['ALT_1/ALT_1'].to_numpy()
    het_nonmissing = (counts['ALT_1/REF'] + counts['ALT_1/ALT_2'] + counts['ALT_1/ALT_1']).to_numpy()
    het_total = (counts['ALT_1/REF'] + counts['ALT_1/ALT_2'] + counts['ALT_1/ALT_1'] + counts['ALT_1/MISSING_ENTRY_{}'.format(sample_name)]).to_numpy()
    all_correct = (counts['REF/REF'] + counts['ALT_1/ALT_1'] + counts['ALT_2/ALT_2']).to_numpy()
    all_total = het_total + (counts['REF/ALT_1'] + counts['REF/ALT_2'] + counts['REF/REF'] + counts['REF/MISSING_ENTRY_{}'.format(sample_name)] +
                             counts['ALT_2/ALT_2'] + counts['ALT_2/ALT_1'] + counts['ALT_2/REF'] + counts['ALT_2/MISSING_ENTRY_{}'.format(sample_name)]).to_numpy()
    return(np.vstack([het_correct, het_nonmissing, het_total, all_correct, all_total]))

'''
When testing a single setting
'''
def single_test(truth_vcf, test_vcf, maf_file):
    FNULL = open(os.devnull, 'w')
    res_table = np.zeros([5, len(mafbins)])
    
    for chrom in chrom_list:
        truth = truth_vcf.format(chrom)
        test = test_vcf.format(chrom)
        maf_file = maf_file.format(chrom)
        maf = {int(pos): af if af <= 0.5 else 1 - af for pos, af in np.loadtxt(maf_file, usecols = (1, 2))}
        subprocess.call(r"java -Xmx4g -jar /storage/software/snpEff/SnpSift.jar concordance -v {} {} > snpsift.out".format(truth, test), shell = True, stderr = FNULL)
        res = pd.read_csv('snpsift.out', sep = "\t", index_col = 1, comment = '#', usecols = np.arange(30))
        chrom_res = get_stats(res, maf)
        print('chr{}:'.format(chrom))
        print(chrom_res)
        res_table = res_table + get_stats(res, maf)
    
    print("Het accuracy:")
    print(",".join(["{:.2f}".format(f) for f in res_table[0, :] / res_table[1, :]]))
    print("Common variants: {:.2f}".format(sum(res_table[0, 2:]) / sum(res_table[1, 2:])))
    print("Het raw numbers:")
    print(",".join(["{:.0f} ({:.2f}%)".format(c, p) for c, p in zip(res_table[0, :], res_table[0, :] * 100 / res_table[2, :])]))
    print("Common variants: {:.0f} ({:.2f}%)".format(sum(res_table[0, 2:]), sum(res_table[0, 2:] * 100 / sum(res_table[2, 2:]))))
    print("All raw numbers:")
    print(",".join(["{:.0f} ({:.2f}%)".format(c, p) for c, p in zip(res_table[3, :], res_table[3, :] * 100 / res_table[4, :])]))
    print("Common variants:")
    print(",".join(["{:.0f} ({:.2f}%)".format(c, p) for c, p in zip(res_table[3, 2:], res_table[3, 2:] * 100 / res_table[4, 2:])]))

'''
When testing multiple settings, to print out result from each setting in a separate column
'''
def multiple_test(truth_vcf, test_vcf_list, maf_file):

    FNULL = open(os.devnull, 'w')
    res_table = np.zeros([len(test_vcf_list), 5, len(mafbins)])
    
    for chrom in chrom_list:
        truth = truth_vcf.format(chrom)
        maf_list = maf_file.format(chrom)
        maf = {int(pos): af if af <= 0.5 else 1 - af for pos, af in np.loadtxt(maf_list, usecols = (1, 2))}
        
        for i, test_vcf in enumerate(test_vcf_list):
            test = test_vcf.format(chrom)
            subprocess.call(r"java -Xmx4g -jar /storage/software/snpEff/SnpSift.jar concordance -v {} {} > snpsift.out".format(truth, test), shell = True, stderr = FNULL)
            res = pd.read_csv('snpsift.out', sep = "\t", index_col = 1, comment = '#', usecols = np.arange(30))
            chrom_res = get_stats(res, maf)
            print('{}, chr{}:'.format(test, chrom))
            print(chrom_res)
            res_table[i, :, :] = res_table[i, :, :] + get_stats(res, maf)
        print(res_table)
            
    print("Het accuracy:")
    for m in range(res_table.shape[2]):
        print(",".join(["{:.2f}".format(f) for f in res_table[:, 0, m] / res_table[:, 1, m]]))

    print("Common variants: " + ",".join(["{:.2f}".format(f) for f in np.sum(res_table[:, 0, 2:], 1) / np.sum(res_table[:, 1, 2:], 1)]))
          
    print("Het raw numbers:")
    for m in range(res_table.shape[2]):
        print(",".join(["{:.0f} ({:.2f}%)".format(c, p) for c, p in zip(res_table[:, 0, m], res_table[:, 0, m] * 100 / res_table[:, 2, m])]))

    print("Common variants: " + ",".join(["{:.0f} ({:.2f}%)".format(c, p) for c, p in zip(np.sum(res_table[:, 0, 2:], 1), np.sum(res_table[:, 0, 2:], 1) * 100 / np.sum(res_table[:, 2, 2:], 1))]))
    
    print("All raw numbers:")
    for m in range(res_table.shape[2]):
        print(",".join(["{:.0f} ({:.2f}%)".format(c, p) for c, p in zip(res_table[:, 3, m], res_table[:, 3, m] * 100 / res_table[:, 4, m])]))

    print("Common variants: " + ",".join(["{:.0f} ({:.2f}%)".format(c, p) for c, p in zip(np.sum(res_table[:, 3, 2:], 1), np.sum(res_table[:, 3, 2:], 1) * 100 / np.sum(res_table[:, 4, 2:], 1))]))
    
multiple_test(truth_vcf, test_vcf_list, maf_file)
