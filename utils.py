import pandas as pd 
from scipy.stats import fisher_exact
import numpy as np 
import os 
import pymysql
from itertools import combinations

from dotenv import load_dotenv
import time


load_dotenv()



def get_query(k,sites):
    sites = tuple(sites)
    if len(sites) == 1:
        q = f'''SELECT exp_condition , mapped_phosphosite FROM profile_insert_table WHERE mapped_gene = '{k}' AND mapped_phosphosite IN ('{sites[0]}') '''
    else:
        q = f'''SELECT exp_condition , mapped_phosphosite FROM profile_insert_table WHERE mapped_gene = '{k}' AND mapped_phosphosite IN {sites}  '''

    print(q)
    return q

def get_query_d(k,e):

    q = f'''SELECT exp_condition , mapped_phosphosite, expression FROM phospodb_nisar_differential_data WHERE mapped_genesymbol ='{k}'  '''

    return q


def get_d(query):
    connection = pymysql.connect(
        host=DB_HOST,
        port=3306,
        user=DB_USER,
        password=DB_PASSWORD,
        database=DB_NAME,
    )
    cursor = connection.cursor()
    cursor.execute(query)
    result = cursor.fetchall()
    df_result = pd.DataFrame(result, columns=[desc[0] for desc in cursor.description])
    cursor.close()
    connection.close()
    return df_result


def get_combinations(df):
    pairs_list = []

    for i in range(len(df)):
        current_data = df.at[i, 'mapped_phosphosite']
        remaining_data = df['mapped_phosphosite'].iloc[i+1:].tolist()
        
        pairs = [(current_data, item) for item in remaining_data]
        
        pairs_list.extend(pairs)

    result_df = pd.DataFrame(pairs_list, columns=['site1', 'site2'])

    mask = result_df['site1'] != result_df['site2']

    result_df = result_df[mask]

    return result_df


def get_n00(s1,s2,t):
    return len(t - (s1.union(s2)))

def get_n01(s1,s2,t):
    return len(s2 - (s1.intersection(s2)))

def get_n10(s1,s2,t):
    return len(s1 - (s1.intersection(s2)))

def get_n11(s1,s2,t):
    return len(s1.intersection(s2))

def get_n(s1,s2,df, total):

    s1_con = df.at[s1,'exp_condition']
    s2_con = df.at[s2,'exp_condition']

    s1_con = set(s1_con)
    s2_con = set(s2_con)
    total = set(total)

    n_00 = get_n00(s1_con,s2_con,total) 
    n_01 = get_n01(s1_con,s2_con,total)
    n_10 = get_n10(s1_con,s2_con,total)
    n_11 = get_n11(s1_con,s2_con,total)
    data = np.array([[n_00,n_01],[n_10,n_11]])
    _ ,p_value = fisher_exact(data, alternative = 'greater')
    return n_00,n_01,n_10,n_11 , p_value

def generate_matrix(df, CUT_OFF,sites):
    
    total_ex = df["exp_condition"].unique().tolist() 
    print(len(total_ex))
    df = df.groupby('mapped_phosphosite').agg(pd.Series.tolist).reset_index()

    df['count'] = df["exp_condition"].apply(lambda x:len(x))
    df = df.sort_values(by=['count'], ascending=False)


    df = df.loc[df['count'] >= CUT_OFF]

    all_sites = pd.DataFrame(df['mapped_phosphosite'])
    all_sites.reset_index(inplace=True)
    df.set_index('mapped_phosphosite', inplace = True)
    print("========================================================================")
    print(all_sites)
    df_comb = get_combinations(all_sites)
    print(df_comb)
    print("========================================================================")


    df_comb[['n_00','n_01','n_10','n_11','p-Value']] = df_comb.apply(lambda x:get_n(x['site1'],x['site2'], df , total_ex), axis = 1, result_type='expand')

    # df_comb[['n_00','n_01','n_10','n_11','p-Value']] = pd.DataFrame(result.tolist())

    df_comb.drop_duplicates(inplace=True)
    df_comb.sort_values(by=['p-Value'], inplace=True)

    return df_comb


def getCDF(df):
    df['CDF'] = np.arange(1, len(df) + 1) / len(df)
    return df

def get_trypsin(s1,s2, df):
    s1 = int(s1[1:])
    s2 = int(s2[1:])

    peptides = []
    for i,row in df.iterrows():
        if (s1 in row['counts']) and (s2 in row['counts']):
            return [row['Fragments']]
        
        elif s1 in row['counts']:
            peptides.append(row['Fragments'])

        elif s2 in row['counts']:
            peptides.append(row['Fragments'])

        else:
            continue
    
    return peptides
    

def trypsin_digest(df):
    
    df_tr = pd.read_excel(r"C:\Users\Admin\Downloads\trypsin_digestion_PAK1.xlsx")
    main_count = []
    couts = []
    c = 1
    for i, row in df_tr.iterrows():
        for _ in range (0, len(row['Fragments'])):
            couts.append(c)
            c += 1
        main_count.append(couts)
        couts = []
    df_tr["counts"] = main_count

    df['peptides'] = df.apply(lambda x:get_trypsin(x['site1'], x['site2'], df_tr), axis =1)

    df["same_peptide"] = df['peptides'].apply(lambda x:"YES" if len(x) ==1 else "NO")
   
    return df






















