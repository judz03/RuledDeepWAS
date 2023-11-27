# Import libraries
from Bio import SeqIO
from Bio.Seq import MutableSeq, Seq
import pandas as pd
import subprocess

def generate_bed(data_path: str, chromosome: str, res_path: str):
    '''
    Function to process the data from the VCF files from Ensembl Variation.
    Limited to single alternate alleles and SNPs. Check the `data_processing.ipynb` file.

    Variables:
    ----------
    data_path: path where the tsv files extracted by bcftools from the original vcf files are stored

    res_path: path where the produced files will be saved

    chromosome: number of chromosome to process
    '''
    
    data_path = data_path
    res_path = res_path
    chromosome = chromosome
    chr_n_data = pd.read_csv(data_path+'chr{}_data.tsv'.format(chromosome), sep = '\t')
    chr_n_data.rename(columns={'#[1]CHROM':'chr', '[2]POS':'pos', '[3]REF':'ref', '[4]ALT':'alt', '[5]TSA':'tsa', '[6]ID':'id'}, inplace=True)
    
    # Filter the variants that are snps and that have only one alternate allele
    chr_n_snps = chr_n_data.loc[(chr_n_data['tsa']=='SNV') & (chr_n_data['alt'].str.len() == 1)]

    chr_n_snps['start'] = chr_n_snps['pos'].astype(int) - 64
    chr_n_snps['end'] = chr_n_snps['pos'].astype(int) + 63
    chr_n_snps['bed'] = chr_n_snps['chr'].astype(str) + ':' + chr_n_snps['start'].astype(str) + '-' + chr_n_snps['end'].astype(str)
    chr_n_snps['bed'].to_csv('{}/bed_chr{}'.format(res_path, chromosome), sep='\t', index=False, header=False)
    return chr_n_snps


generate
    # This code has to be run directly in the terminal after adding the samtools path to the PATH variable
    exit_code = subprocess.Popen("samtools faidx /mnt/sda1/Databases/Reference Genome/GRCh38p14/Ensembl/Homo_sapiens_GRCh38_dna_primary_assembly.fa -r {}/bed_chr{} -o {}/seq_vcf_chr{}".format(res_path, chromosome, res_path, chromosome), 
                                 shell=True, stdout=subprocess.PIPE).stdout.read()
    if str(exit_code,'utf-8')!='':
        print(str(exit_code,'utf-8'))
        exit(1)
    
    records = list(SeqIO.parse(res_path+'seq_vcf_chr{}'.format(chromosome), "fasta"))
    ref_seqs = [sequence[1].seq for sequence in enumerate(records)]

    alt_alleles = list(chr_n_snps['alt'])
    alt_seqs = []
    i = 0
    while i < len(ref_seqs):
        mutable = MutableSeq(ref_seqs[i])
        mutable[64] = alt_alleles[i]
        alt_seqs.append(str(mutable))
        i+=1
    
    ref_seqs = list(map(str, ref_seqs)) # Turn Seq objects into normal Python Strings

    chr_n_df = pd.DataFrame(data = chr_n_snps, copy=True)
    chr_n_df['ref_seq'] = ref_seqs
    chr_n_df['alt_seq'] = alt_seqs

    return chr_n_df
