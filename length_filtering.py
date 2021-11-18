import os
import pathlib
import shutil
import sys
import tempfile

import pandas as pd
from Bio import SeqIO
from qiime2 import Artifact, Metadata
from qiime2.plugins.feature_table.methods import filter_features, filter_seqs


    #Length filter of rep_seqs
def length_filter(representative_sequences, feature_table, length_to_filter):
    #Pull fasta from rep_seqs to use as metadata
    def extract_fasta(file, dest):
        with tempfile.TemporaryDirectory() as temp:
            file.export_data(temp)
            temp_pathlib = pathlib.Path(temp)
            for item in temp_pathlib.iterdir():
                if item.suffix == '.fasta':
                    shutil.copy(item, dest)

    os.system('mkdir fastas')
    extract_fasta(representative_sequences, 'fastas')

    #Parse fasta to get IDs and lengths, then send to df and filter IDs by length
    with open('fastas/dna-sequences.fasta') as fasta_file:
        featureids = []
        lengths = []
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):
            featureids.append(seq_record.id)
            lengths.append(len(seq_record.seq))

    rep_seqs_meta = pd.DataFrame({"FeatureID":featureids, "Length":lengths})

    features_to_exclude = rep_seqs_meta[rep_seqs_meta['Length'] < length_to_filter]

    features_to_exclude['FeatureID'].to_csv('Features-to-exclude.csv', index=False)


    exclude = Metadata.load("Features-to-exclude.csv")

    #Filter rep-seqs based on seqs_to_exclude
    rep_seqs_filt = filter_seqs(representative_sequences,
                               metadata = exclude,
                               exclude_ids = True)

    #Filter table based on seqs_to_exclude
    table_filt = filter_features(feature_table,
                                metadata = exclude,
                                exclude_ids = True)

    return rep_seqs_filt, table_filt

if __name__ == '__main__':
    main()

