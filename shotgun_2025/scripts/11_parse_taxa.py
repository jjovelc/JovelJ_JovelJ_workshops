import pandas as pd
import sys

# Load the data
input_file = sys.argv[1]
df = pd.read_csv(input_file, sep='\t', index_col='ID')

# Function to extract and save subtables for each taxonomic level
def save_subtables(df, base_filepath):
    taxa_split = df.index.str.split('|')
    max_depth = max(len(taxa) for taxa in taxa_split)
    
    for depth in range(max_depth):
        taxa_level = ['|'.join(taxa[:depth + 1]) if len(taxa) > depth else '' for taxa in taxa_split]
        sub_df = df.copy()
        sub_df.index = taxa_level
        sub_df = sub_df[~sub_df.index.duplicated(keep='first')]  # Remove duplicates
        sub_df.reset_index(inplace=True)
        sub_df.rename(columns={'index': 'Taxa'}, inplace=True)
        
        # Filter out rows where 'Taxa' is an empty string
        sub_df = sub_df[sub_df['Taxa'] != '']
        
        taxonomic_rank = f'level_{depth + 1}'
        output_filename = f'{base_filepath}_{taxonomic_rank}.tsv'
        sub_df.to_csv(output_filename, sep='\t', index=False)
        print(f'Saved: {output_filename}')

# Extract the base filepath without extension to use in output filenames
base_filepath = input_file.rsplit('.', 1)[0]

save_subtables(df, base_filepath)
