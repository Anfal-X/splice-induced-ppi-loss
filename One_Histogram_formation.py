import pandas as pd
import matplotlib.pyplot as plt

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

# Create one figure for all overlaid curves
fig, ax = plt.subplots(figsize=(12, 7))

# Loop through each model file
for file, model in zip(files, model_names):

    # Load the file
    table = pd.read_csv(file)

    results = []

    # Only negative isoforms — those that lost the interaction
    for index, row in table[table['Interaction_Outcome'] == 'negative'].iterrows():

        spliced_positions = set(eval(str(row['Spliced_Positions'])))
        all_scores = eval(str(row['Full_PiPENN_Scores']))

        if len(spliced_positions) == 0:
            continue

        # Positions predicted as interface
        predicted_interface = set([i for i, score in enumerate(all_scores) if score > 0.5])

        # Total predicted interface size
        interface_size = len(predicted_interface)

        # Number of interface residues spliced out
        interface_lost = len(spliced_positions & predicted_interface)

        results.append({
            'Isoform_ID': row['Isoform_ID'],
            'Predicted_Interface_Size': interface_size,
            'Interface_Lost': interface_lost
        })

    results_df = pd.DataFrame(results)
    print(f'{model}: {len(results_df)} negative isoforms analysed')

    # Save individual results CSV for each model
    results_df.to_csv(f'interface_lost_results_{model}.csv', index=False)

    # Plot as a density curve overlaid on the same figure
    
    results_df['Interface_Lost'].plot.kde(ax=ax, label=model, linewidth=2)

ax.set_title('Interface Residues Lost by Splicing Across All Models\n(Negative Isoforms Only)', fontsize=14)
ax.set_xlabel('Number of interface residues spliced out', fontsize=12)
ax.set_ylabel('Density (n=207 per model)', fontsize=12)
ax.legend(title='Model', fontsize=9)
ax.set_xlim(left=0)

plt.tight_layout()
plt.savefig('interface_lost_all_models_overlaid.png', dpi=300)
plt.show()

print('Done!')