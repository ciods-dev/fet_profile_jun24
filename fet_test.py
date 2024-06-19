import pandas as pd
from utils import get_query, get_d, generate_matrix

from myplot import get_plot

def getfet(kianase, CUT_OFF,sites):

    q = get_query(kianase,sites)
    df = get_d(q)
    df = generate_matrix(df , CUT_OFF ,sites)

    return df

df = kinase_data = pd.read_excel("GetTOP10.xlsx",)

kinase  = df["HGNC Name"].unique()

final_df = pd.DataFrame()
single_df = pd.DataFrame()

for k in kinase:
    sites = df.loc[ df["HGNC Name"] == k , "mapped_phosphosite"].tolist()
    if len(sites) > 1:
        fetdf = getfet(k,0,sites)
        fetdf["kinase"] = k
        final_df = pd.concat([final_df,fetdf])
    elif len(sites) == 1:
        s_sites = df.loc[ df["HGNC Name"] == k ]

        single_df = pd.concat([single_df,s_sites])


final_df.to_excel("Final_result_sheet_1.xlsx" )
single_df.to_excel("single_df.xlsx" )