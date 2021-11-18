import pathlib
import shutil
import tempfile
import os
import numpy as np
import biom
import pandas as pd
from qiime2 import Artifact, Metadata


def per_sample_filter(filtering_integer, feature_table):


    #Importing from a biom table
    try:
        working_table = biom.load_table(feature_table)
    except:
        print('Could not import from biom table, trying from qza')
        try:
    #Importing from imported qza
            working_table = feature_table.view(biom.Table)
        except:
            try:
    #Importing from converted qza
                working_table = (feature_table.filtered_table).view(biom.Table)
            except:
                print('Could not import from qiime artifact, check to ensure input file is appropriate')



    #Convert table to df
    biom_df = working_table.to_dataframe()
    print('File loaded! Filtering at '+str(filtering_integer)+' occurences per feature')

    #Filter out any feature from a sample if it has less than X reads,
    biom_df = biom_df._get_numeric_data()
    biom_dense = biom_df.sparse.to_dense()

    biom_dense[biom_dense < filtering_integer] = 0


    #Translate the filtered table back into a qza
    biom_dense = biom_dense.transpose()
    freq_filt_table = Artifact.import_data("FeatureTable[Frequency]", biom_dense)

    #Export the new table as a qza for further use in Qiime2
    freq_filt_table.save('freq_filt_table.qza')

    return freq_filt_table


if __name__ == '__main__':
    main()
