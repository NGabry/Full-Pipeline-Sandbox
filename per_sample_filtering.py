import pathlib
import shutil
import tempfile
import os
import numpy as np
import biom
import pandas as pd
from qiime2.plugins.feature_table.visualizers import summarize as summarize_table
from qiime2 import Artifact, Metadata


def per_sample_filter(filtering_integer, feature_table):

    #Importing from a biom table
    try:
        table = biom.load_table(feature_table)
    except:
        pass

    #Or importing from a qza file and translating to a biom table
    try:
        qza_table = Artifact.load(feature_table)
        table = qza_table.view(biom.Table)
    except:
        pass

    #Convert table to df
    biom_df = table.to_dataframe()
    biom_df

    #Filter out any feature from a sample if it has less than X reads,
    biom_df = biom_df._get_numeric_data()
    biom_dense = biom_df.sparse.to_dense()

    biom_dense[biom_dense < filtering_integer] = 0


    #Translate the filtered table back into a qza
    biom_dense = biom_dense.transpose()
    freq_filt_table = Artifact.import_data("FeatureTable[Frequency]", biom_dense)

    #Export the new table as a qza for further use in Qiime2
    freq_filt_table.save('freq_filt_table.qza')


if __name__ == '__main__':
    main()
