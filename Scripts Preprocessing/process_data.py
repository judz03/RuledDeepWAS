# Import libraries
from Bio import SeqIO
from Bio.Seq import MutableSeq, Seq
import pandas as pd
import subprocess

def generate_bed(chr_n_snps: pd.DataFrame, len_sequence: int = 64, chromosome: str=None, res_path: str = None, gen_bed:bool = False):
    '''
    Function to process the data from the VCF files from Ensembl Variation.
    Limited to single alternate alleles and SNPs. Check the `data_processing.ipynb` file for additional reference

    Variables:
    ----------
    chr_n_snps: dataframe containing the variant registers (only SNVs)

    res_path: path where the produced files will be saved (required if gen_bed = True)

    chromosome: number of chromosome to process (required if gen_bed = True)

    gen_bed: A boolean value to indicate if the `bed` column will be exported into a csv file
    '''
    
    # Create the list of alternative alleles
    chr_n_snps["alt_list"] = chr_n_snps["alt"].str.split(pat=",")
    chr_n_snps['start'] = chr_n_snps['pos'].astype(int) - (len_sequence-1)
    chr_n_snps['end'] = chr_n_snps['pos'].astype(int) + len_sequence
    chr_n_snps['bed'] = chr_n_snps['chr'].astype(str) + ':' + chr_n_snps['start'].astype(str) + '-' + chr_n_snps['end'].astype(str)
    if gen_bed == True: chr_n_snps['bed'].to_csv('{}/bed_chr{}'.format(res_path, chromosome), sep='\t', index=False, header=False)
    return chr_n_snps





def generate_sequences(variant_df: pd.DataFrame, chromosome: str, seq_path: str, bed_path: str=None, ref_genome_path:str=None,
                       generate_fasta: bool = False):
    '''
    This function generates the reference and alternative sequences based on the information in variant_df and the reference genome. Calls samtools to do this
    later task, and then saves the generated sequences in `gwas_associated_seq_path`.

    Parameters:
    ----------------
    variant_df: a DataFrame containing the genomic variants information. These are the result of extracting the information from the original VCF files with bcftools.

    ref_genome_path: ()

    seq_path

    bed_path

    chromosome

    generate_fasta
    
    '''
    
    if generate_fasta == True:
        # This code has to be run directly in the terminal after adding the samtools path to the PATH variable
        # These lines were extracted from the DeepPerVar paper
    #-----------------------------------------------------------------------------------------------------#
        exit_code = subprocess.Popen("samtools faidx {} -r {}/bed_chr{} -o {}/ref_seq_chr{}".format(ref_genome_path, bed_path,
                                                                                                    chromosome, seq_path, chromosome),
                                    shell=True, stdout=subprocess.PIPE).stdout.read()
        if str(exit_code,'utf-8')!='':
            print(str(exit_code,'utf-8'))
            exit(1)

    records = list(SeqIO.parse(seq_path+'/ref_seq_chr{}'.format(chromosome), "fasta"))
    ref_seqs = [sequence[1].seq for sequence in enumerate(records)]
    ref_seqs = list(map(str, ref_seqs)) # Turn Seq objects into normal Python Strings and save them in a list
    variant_df['ref_seq'] = ref_seqs
    #-----------------------------------------------------------------------------------------------------#
    # Generate the subsequences for variants, including the multiallelic ones
    alt_seq = []
    for idx, variant in enumerate(variant_df['alt_list']):
        tmp_alt_seq = [] # Clear the contents of this list each time the for loop goes to a new register
        for allele in variant:
            tmp_seq = MutableSeq(variant_df['ref_seq'].iloc[idx])
            tmp_seq[63] = allele
            tmp_seq = str(tmp_seq)
            tmp_alt_seq.append(tmp_seq)
        alt_seq.append(tmp_alt_seq)

    variant_df['alt_seq'] = alt_seq
    return variant_df