import pandas as pd 
import matplotlib.pyplot as plt
#How much of the protein interface is spliced out that interaction is changed from positive to negative
File=pd.read_csv('combined_web_dnet_input_H.csv')# this file is contains combined content from longest and reference 1 final tables
results=[] # empty file to store results 
# Filters the rows where interaction outcome is negative
for index, row in File[File['Interaction_Outcome']=='negative'].iterrows():
    s_p=set(eval(str(row['Spliced_Positions'])))# reads spliced positions for the specific isoform; set() used to make overlap calculation easier
    pipenn_scores=eval(str(row['Full_PiPENN_Scores'])) # reads pipenn scores 
    if len(s_p)==0:# skips isoforms who do not have any spliced postions
        continue
    # provides all residues position which are predicted as interface one having pipenn scores greater than 0.5
    all_residues_predicted_interface = set([i for i, score in enumerate(pipenn_scores) if score > 0.5])
    no_of_residues_as_interface=len(all_residues_predicted_interface)# total number of residues predicted as interface for the specific isoform
    # Number of lost residues which were predicted as interface;intersection between spliced positions and interaction interface
    no_of_residues_spliced_out = len(s_p & all_residues_predicted_interface)
    # calculates percentage of residues which were part of interface and were lost during splicing
    if no_of_residues_as_interface > 0:
        percent_lost = (no_of_residues_spliced_out / no_of_residues_as_interface) * 100
    else:
        percent_lost = 0
    results.append({
        'Isoform_ID':row['Isoform_ID'],
        'Percentage_of_interface_lost': percent_lost
    })
results_df=pd.DataFrame(results) # converts list into table
results_df.to_csv('Percentage_of_interface_lost_per_isoform.csv',index=False)

# Building the histogram
plt.figure(figsize=(8, 5)) # creates a figure
plt.hist(results_df['Percentage_of_interface_lost'], bins=20, color='blue', edgecolor='black')
# adds axis and title
plt.title('Percentage of Predicted Interface Lost by Splicing')
plt.xlabel('% of interface spliced out')
plt.ylabel('Number of isoforms')
plt.tight_layout() # nothing overlaps or get cutoff
plt.savefig('ensemble_mixed_pairs_interface_lost_histogram.png', dpi=300)
plt.show()
print('Done!')