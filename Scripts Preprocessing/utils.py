'''
This is a library that contains multipleuseful functions for the preprocessing of the databases used
in the project.

Functions:
----------
* replace_seq_info
* export_variants_df
'''
# Import the necessary modules and libraries
from Bio import SeqIO
# To be able to modify sequences
from Bio.Seq import MutableSeq, Seq
# Import SeqRecord
from Bio.SeqRecord import SeqRecord
# Import pandas
import pandas as pd

# %% Replace the variants information in the reference sequence
def replace_seq_info(reference_sequence:MutableSeq, variant_info:pd.Series, context:int=128)->SeqRecord:
    # We will have to import Bio.SeqIO and Bio.Seq.MutableSeq
    '''
    This function replaces the information in the reference sequence with the specifications
    of the variant databases.
    **** LIMITED TO SNVS ****

    Parameters:
    -----------
    reference_sequence: A MutableSeq object loaded trough Biopython which contains the reference
    sequence of a chromosome or genome.

    variant_info: Relevant characteristics about the variant which include:

        * variant_seq: The information of the variant allele which contains the nucleotide(s) to 
        replace in the reference sequence.
        
        * start: The coordenate where the variation in the reference sequence begins.

        * end: The coordenate where the variation in the reference sequence ends.

    context: The number of bp after and before the start and end, respectively, to be considered
    to produce a subsequence with context surrounding the variant.
    '''
    mutable_sequence = reference_sequence
    variant_info = variant_info
    start = int(variant_info['start'])-1
    end = int(variant_info['end']) # Se quito un -1
    # Insert mutation
    mutable_sequence[start] = variant_info['Variant_seq']
    # Extract subsequence with mutation and context
    variant_sequence = mutable_sequence[start-context:end+context] #Se quito un +1 despues de context
    variant_sequence = variant_sequence.upper()
    # Convert into unmutable sequence
    variant_sequence = Seq(variant_sequence)
    # Turn into a SeqRecord object and associate relevant information
    variant_sequence = SeqRecord(variant_sequence)
    variant_sequence.id = variant_info['ID']
    variant_sequence.name = 'Sequence containing the variant with ID: ' + variant_info['ID']
    variant_sequence.dbxrefs = [variant_info['Dbxref']]
    variant_sequence.annotations['Chromosome'] = variant_info['seqid']
    variant_sequence.annotations['Clinical_significance'] = variant_info['clinical_significance']
    variant_sequence.annotations['Start of mutation'] = variant_info['start']
    variant_sequence.annotations['End of mutation'] = variant_info['end']
    return variant_sequence


def export_variants_df(mutable_sequence: MutableSeq|SeqRecord|Seq, variants_info:pd.DataFrame)-> pd.DataFrame:
    '''
    This function exports the modified sequences with all the relevant information to a DataFrame.

    Parameters:
    -----------
    mutable_sequence: the nucleotide sequence to be modified, can be a MutableSeq, SeqRecord or Seq
    
    '''
    if type(mutable_sequence) == SeqRecord:
        mutable_sequence = MutableSeq(mutable_sequence.seq)
    elif type(mutable_sequence) == Seq:
        mutable_sequence = MutableSeq(mutable_sequence)

    variants_info = variants_info
    mutable_sequence = mutable_sequence
    variant_sequences_data = {'sequence':[],
                          'id':[],
                          'chromosome':[],
                          'start':[],
                          'end':[],
                          'Dbxref':[],
                          'clinical_significance':[]}

    for index in range(len(variants_info)):
        variant_sequence_info = replace_seq_info(mutable_sequence,variants_info.iloc[index])
        variant_sequences_data['sequence'].append(variant_sequence_info.seq)
        variant_sequences_data['id'].append(variant_sequence_info.id)
        variant_sequences_data['chromosome'].append(variant_sequence_info.annotations['Chromosome'])
        variant_sequences_data['start'].append(variant_sequence_info.annotations['Start of mutation'])
        variant_sequences_data['end'].append(variant_sequence_info.annotations['End of mutation'])
        variant_sequences_data['Dbxref'].append(variant_sequence_info.dbxrefs)
        variant_sequences_data['clinical_significance'].append(variant_sequence_info.annotations['Clinical_significance'])

    df = pd.DataFrame(variant_sequences_data)
    return df