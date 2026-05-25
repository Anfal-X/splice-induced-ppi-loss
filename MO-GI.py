import pandas as pd
#File='combined_ensemble_input_H.csv'
# All 8 model files
files = [
    'combined_ensemble_input_H.csv',
    'combined_rnet_input_H.csv',
    'combined_cnn_rnn_input_H.csv',
    'combined_dnet_input_H.csv',
    'combined_ann_input_H.csv',
    'combined_unet_input_H.csv',
    'combined_rnn_input_H.csv',
    'combined_web_dnet_input_H.csv'
]

# Clean model names for the legend
model_names = [
    'Ensemble',
    'RNET',
    'CNN-RNN',
    'DNET',
    'ANN',
    'UNET',
    'RNN',
    'Web-DNET'
]
for file, model in zip(files, model_names):
    output=pd.read_csv(file)
    # groups of unique gene and interactor combination are made
    pairs = output.groupby(['Gene_Symbol', 'Interactor_Symbol'])
    different_outcome_pairs=[]
    for (gene, interactor), group in pairs:# loops over every gene-interactor pairs
        outcomes = group['Interaction_Outcome'].unique()# takes out unique interaction outcome in this group
        if 'negative' in outcomes and 'positive' in outcomes:# checks wether this group has one negative and one positive outcome
            different_outcome_pairs.append(group)
            print(f'Mixed pair found: {gene} - {interactor} ({len(group)} isoforms)')

    if different_outcome_pairs:
        different_outcomes_df = pd.concat(different_outcome_pairs)
        #output_name = File.replace('.csv', '_mixed_pairs.csv')
        different_outcomes_df.to_csv(f'Triples_{model}.csv', index=False)
        print(f'\nTotal isoforms in  pairs: {len(different_outcomes_df)}')
        print(f'Total unique outcome pairs found: {len(different_outcome_pairs)}')
        print(different_outcomes_df[['Isoform_ID', 'Gene_Symbol', 'Interactor_Symbol', 'Interaction_Outcome']])
    else:
        print('No mixed outcome pairs found')       