import pandas as pd
from scipy import stats
File='Triples_Ensemble.csv'
table = pd.read_csv(File)
print('Rows loaded:', len(table))
print(table['Interaction_Outcome'].value_counts())
# Means scores contains numbers
table['Mean_Scores']=pd.to_numeric(table['Mean_Scores'],errors='coerce')# convert anything into number which is saved as text
# splitting of the table into two groups containing positive and negative scores
negative_s=table[table['Interaction_Outcome']=='negative']['Mean_Scores'].dropna()# removes any missing values
positive_s=table[table['Interaction_Outcome']=='positive']['Mean_Scores'].dropna()
print('Negative group:',len(negative_s))
print('Postive group:',len(positive_s))
print('Mean negative score:',round(negative_s.mean(),4))
print('Mean posituve score:',round(positive_s.mean(),4))
print('--- Wilcoxon rank-sum test ---')
statistic, pvalue = stats.mannwhitneyu(negative_s, positive_s, alternative='two-sided')
print('Statistic:', round(statistic, 4))
print('P-value:', pvalue)


print('--- Jaccard index ---')
jaccard_results = []
for index,row in table.iterrows():# loops over every row
    isoform_id = row['Isoform_ID']
    outcome = row['Interaction_Outcome']
    #converts back into actual python list as when pandas reads a list in csv it stores it as string
    spliced_positions = eval(str(row['Spliced_Positions']))
    All_pipenn_scores = eval(str(row['Full_PiPENN_Scores']))
    # gives set of all positions predicted as interface in the protein
    predicted_interface_positions = set([i for i, score in enumerate(All_pipenn_scores) if score > 0.5])
    spliced_set = set(spliced_positions)#converts list into sets
    intersection = len(spliced_set & predicted_interface_positions)# includes positions that are in both sets
    union = len(spliced_set |predicted_interface_positions)# include positions that are in one of the each set
    if union == 0:
        continue
    jaccard = intersection / union
    jaccard_results.append({
        'Isoform_ID': isoform_id,
        'Interaction_Outcome': outcome,
        'Jaccard_Score': jaccard
    })
jaccard_df = pd.DataFrame(jaccard_results)
print('Jaccard rows:', len(jaccard_df))
print(jaccard_df['Interaction_Outcome'].value_counts())
print('Mean Jaccard negative:', round(jaccard_df[jaccard_df['Interaction_Outcome'] == 'negative']['Jaccard_Score'].mean(), 4))
print('Mean Jaccard positive:', round(jaccard_df[jaccard_df['Interaction_Outcome'] == 'positive']['Jaccard_Score'].mean(), 4))
neg_jaccard = jaccard_df[jaccard_df['Interaction_Outcome'] == 'negative']['Jaccard_Score']
pos_jaccard = jaccard_df[jaccard_df['Interaction_Outcome'] == 'positive']['Jaccard_Score']
jaccard_stat, jaccard_pvalue = stats.mannwhitneyu(neg_jaccard, pos_jaccard, alternative='two-sided')
print('Jaccard Wilcoxon statistic:', round(jaccard_stat, 4))
print('Jaccard Wilcoxon p-value:', jaccard_pvalue)


print('--- Chi-square test ---')
spliced_and_predicted_interface = 0
spliced_and_not_predicted_interface = 0
not_spliced_and_predicted_interface = 0
not_spliced_and_not_predicted_interface = 0
for index, row in table.iterrows():# goes through every row in the final table
    spliced_positions = set(eval(str(row['Spliced_Positions'])))# reads the spliced positions from that row
    all_scores = eval(str(row['Full_PiPENN_Scores'])) #list of pipenn scores
    all_positions = set(range(len(all_scores)))# sequence of numbers containing every residue position in protein
    interface_positions = set([i for i, score in enumerate(all_scores) if score > 0.5])# marking positions which are predicted as interface
    for position in all_positions:# for every residue it is being checked whether it is predicted interface or spliced out residue and marks True/False
        is_spliced = position in spliced_positions
        is_interface = position in interface_positions
        if is_spliced and is_interface:
            spliced_and_predicted_interface += 1
        elif is_spliced and not is_interface:
            spliced_and_not_predicted_interface += 1
        elif not is_spliced and is_interface:
            not_spliced_and_predicted_interface += 1
        else:
            not_spliced_and_not_predicted_interface += 1
        contingency_table = [[spliced_and_predicted_interface, spliced_and_not_predicted_interface],
                     [not_spliced_and_predicted_interface, not_spliced_and_not_predicted_interface]]

chi2, pvalue_chi2, dof, expected = stats.chi2_contingency(contingency_table)
print('Contingency table:')
print('                        Predicted interface    Not predicted interface')
print('Spliced-out:           ', spliced_and_predicted_interface, '                      ', spliced_and_not_predicted_interface)
print('Not spliced-out:       ', not_spliced_and_predicted_interface, '                      ', not_spliced_and_not_predicted_interface)
print('Chi-square statistic:', round(chi2, 4))
print('Degrees of freedom:', dof)
print('P-value:', pvalue_chi2)

jaccard_df.to_csv(File + '_jaccard_scores.csv', index=False)
print('Individual Jaccard scores saved as:', File + '_jaccard_scores.csv')

results_summary = pd.DataFrame([
    {'Test': 'Wilcoxon', 'Comparison': 'Mean PiPENN scores at spliced positions', 'Statistic': round(statistic, 4), 'P_value': pvalue},
    {'Test': 'Jaccard', 'Comparison': 'Mean Jaccard negative vs positive', 'Statistic_negative': round(neg_jaccard.mean(), 4), 'Statistic_positive': round(pos_jaccard.mean(), 4), 'P_value': jaccard_pvalue},
    {'Test': 'Chi-square', 'Comparison': 'Spliced-out vs predicted interface residue level', 'Statistic': round(chi2, 4), 'P_value': pvalue_chi2}
])
contingency_df = pd.DataFrame(
   [[spliced_and_predicted_interface, spliced_and_not_predicted_interface],
    [not_spliced_and_predicted_interface, not_spliced_and_not_predicted_interface]],
    index=['Spliced_out', 'Not_spliced_out'],
    columns=['Predicted_interface', 'Not_predicted_interface']
)
contingency_df.to_csv(File + '_contingency_table.csv', index=True)
print('Contingency table saved as:', File + '_contingency_table.csv')

results_summary.to_csv(File + '_statistical_results.csv', index=False)
#print('Results saved as:', File + '_statistical_results.csv')
# Save p-values to a small results file
import os
result_row = pd.DataFrame([{
    'File': File,
    'Wilcoxon_p': pvalue,
    'Jaccard_p': jaccard_pvalue,
    'Chi2_p': pvalue_chi2
}])
#Check if the output file already exists
#If yes — append to it. If no — create it with a header
output_file = 'all_pvalues.csv'
if os.path.exists(output_file):
    result_row.to_csv(output_file, mode='a', header=False, index=False)
    #mode='a' means append — adds a new row without overwriting
else:
    result_row.to_csv(output_file, mode='w', header=True, index=False)
    #mode='w' means write fresh — creates the file with column names
print(f'P-values appended to {output_file}')