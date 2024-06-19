import pandas as pd 

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


df = pd.DataFrame([['T36'],['S75']], columns= ['mapped_phosphosite']) 

result_df = get_combinations(df)

print(result_df)