import os.path
import glob as glob
import pandas as pd
import natsort as natsort

def find_sample(fname):
    return(os.path.basename(fname))

list_folders = natsort.natsorted(glob.glob("./ALIGNED/*mESC*"))
ref_path = list_folders[0]
df_ref = df = pd.read_table(ref_path+"/ReadsPerGene.out.tab", header=3)
df_ref = df.set_axis(['gene_id', 'R1', 'R2', 'R12'], axis=1)
df_merged = df_ref[['gene_id', 'R1']]
df_merged = df_merged.rename(columns={"R1": find_sample(ref_path)})

for folder in list_folders[1:]:
    sample_name = find_sample(folder)
    print(sample_name)
    path_raw_count = folder+"/ReadsPerGene.out.tab"
    df_tmp = pd.read_table(path_raw_count, header=3)
    df_tmp = df_tmp.set_axis(['gene_id', 'R1', 'R2', 'R12'], axis=1)
    df_tmp_parsed = df_tmp[['gene_id', 'R1']]
    df_tmp_parsed = df_tmp_parsed.rename(columns={"R1": sample_name})
    df_merged = pd.merge(df_merged, df_tmp_parsed, on='gene_id', how='outer')

df_merged.to_csv("./count_matrix_mESC.csv", index=False)
