import pathlib
import sys
import shutil
import tempfile
import os
from qiime2 import Artifact, Metadata
import pandas as pd
from qiime2.plugins import cutadapt as cut
from qiime2.plugins.demux.visualizers import summarize
from qiime2.plugins.feature_table.visualizers import summarize as summarize_table
from qiime2.plugins.dada2.methods import denoise_single
from qiime2.plugins.feature_table.methods import filter_seqs
from qiime2.plugins.feature_table.methods import filter_features
from qiime2.plugins.feature_table.visualizers import tabulate_seqs


def calculate_IJR(seqs, metadata, meta_df):
    def extract_tsvs(viz, dest):
        with tempfile.TemporaryDirectory() as temp:
            viz.export_data(temp)
            temp_pathlib = pathlib.Path(temp)
            for file in temp_pathlib.iterdir():
                if file.suffix == '.tsv':
                    shutil.copy(file, dest)
    demux_summary = summarize(seqs)
    demux_viz     = demux_summary.visualization
    os.system('mkdir tsvs')
    extract_tsvs(demux_summary.visualization, 'tsvs')
    per_sample_df = pd.read_csv("./tsvs/per-sample-fastq-counts.tsv", sep="\t")
    df            = per_sample_df
    df['Fwd']     = df['sample ID'].str[1:3]
    df['Rev']     = df['sample ID'].str[4:6]
    df            = df.sort_values(['sample ID'])
    df.reset_index(inplace=True)
    df['True_False'] = meta_df['True_False']
    F1            = df[df['Fwd'] == '01']
    R11           = df[df['Rev'] == '11']
    F1_False      = F1[F1["True_False"] == 1]
    F1_jump_rate  = (F1_False['forward sequence count'].sum())/(F1['forward sequence count'].sum())
    R11_False     = R11[R11["True_False"] == 1]
    R11_jump_rate = (R11_False['forward sequence count'].sum())/(R11['forward sequence count'].sum())
    jump_rate     = (F1_jump_rate + R11_jump_rate) / 2
    Indexes       = df['sample ID'].tolist()
    Seq_counts    = df['forward sequence count'].tolist()
    Fwds          = df['Fwd'].tolist()
    Fwds          = list(set(Fwds))
    Fwds.sort()
    Revs          = df['Rev'].tolist()
    Revs          = list(set(Revs))
    Revs.sort()
    indexes_exp_false = []
    total_read_count  = df['forward sequence count'].sum()

    count = 0

    for i in Fwds:
        temp_df = df[df['Fwd'] == i]
        i_seqs = temp_df['forward sequence count'].sum()
        for j in Revs:
            temp_df2 = df[df['Rev'] == j]
            j_seqs = (temp_df2['forward sequence count'].sum()) - Seq_counts[count]
            f_to_r = i_seqs * (j_seqs / total_read_count) * jump_rate
            r_to_f = j_seqs * (i_seqs / total_read_count) * jump_rate
            indiv_tag_jump_count = f_to_r + r_to_f
            to_add = (Indexes[count], indiv_tag_jump_count)
            indexes_exp_false.append(to_add)
            count += 1



    out_df = pd.DataFrame(indexes_exp_false, columns = ['SampleIndex', 'Expected False Reads'])
    out_df.to_csv('Expected_False_Read_Per_Index.csv', index=False)
    max_IJR = out_df['Expected False Reads'].max()
    print('Maximum number of false reads expected in a single sample:', max_IJR)
    return(max_IJR)

def recalculate_IJR(rep_seqs, table, meta_df):

    def extract_csvs(viz, dest):
        with tempfile.TemporaryDirectory() as temp:
            viz.export_data(temp)
            temp_pathlib = pathlib.Path(temp)
            for file in temp_pathlib.iterdir():
                if file.suffix == '.csv':
                    shutil.copy(file, dest)

    table_sum = summarize_table(table)
    table_viz = table_sum.visualization
    os.system('mkdir csvs')
    extract_csvs(table_viz, 'csvs')
    per_sample_df = pd.read_csv("./csvs/sample-frequency-detail.csv", sep=",")
    df            = per_sample_df
    df            = df.rename(columns={'Unnamed: 0': 'sample ID', '0': 'forward sequence count'})
    meta_df2      = meta_df.rename(columns={'#SampleID': 'sample ID'})
    df2           = meta_df2.join(df.set_index('sample ID'), on='sample ID')
    df2           = df2.fillna(0)
    df            = df2
    df['Fwd']     = df['sample ID'].str[1:3]
    df['Rev']     = df['sample ID'].str[4:6]
    df            = df.sort_values(['sample ID'])
    df.reset_index(inplace=True)
    df['True_False'] = meta_df['True_False']
    F1            = df[df['Fwd'] == '01']
    R11           = df[df['Rev'] == '11']
    F1_False      = F1[F1["True_False"] == 1]
    F1_jump_rate  = (F1_False['forward sequence count'].sum())/(F1['forward sequence count'].sum())
    R11_False     = R11[R11["True_False"] == 1]
    R11_jump_rate = (R11_False['forward sequence count'].sum())/(R11['forward sequence count'].sum())
    jump_rate     = (F1_jump_rate + R11_jump_rate) / 2
    Indexes       = df['sample ID'].tolist()
    Seq_counts    = df['forward sequence count'].tolist()
    Fwds          = df['Fwd'].tolist()
    Fwds          = list(set(Fwds))
    Fwds.sort()
    Revs          = df['Rev'].tolist()
    Revs          = list(set(Revs))
    Revs.sort()
    indexes_exp_false = []
    total_read_count  = df['forward sequence count'].sum()

    count = 0

    for i in Fwds:
        temp_df = df[df['Fwd'] == i]
        i_seqs = temp_df['forward sequence count'].sum()
        for j in Revs:
            temp_df2 = df[df['Rev'] == j]
            j_seqs = (temp_df2['forward sequence count'].sum()) - Seq_counts[count]
            f_to_r = i_seqs * (j_seqs / total_read_count) * jump_rate
            r_to_f = j_seqs * (i_seqs / total_read_count) * jump_rate
            indiv_tag_jump_count = f_to_r + r_to_f
            to_add = (Indexes[count], indiv_tag_jump_count)
            indexes_exp_false.append(to_add)
            count += 1



    out_df = pd.DataFrame(indexes_exp_false, columns = ['SampleIndex', 'Expected False Reads'])
    out_df.to_csv('Expected_False_Read_Per_Index_Recalc.csv', index=False)
    max_IJR = out_df['Expected False Reads'].max()
    print('Maximum number of false reads expected in a single sample:', max_IJR)
    return(max_IJR)

if __name__ == '__main__':
    main()
