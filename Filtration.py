import pandas as pd
#loading result file of pipenn-emb
protein_predictions=pd.read_csv('preds-ensemble-all.csv')
print('Total_proteins:',len(protein_predictions))
#loading interactor table 2B from yang excel sheet 11 
interactor_panel=pd.read_excel('Fixed_interactor_table.xlsx')
print('interactor_panel rows:',len(interactor_panel))
#Filter 1 : removing of rows where interaction found is N/A
F1=interactor_panel[interactor_panel['Interaction_Found'].isin(['positive','negative'])]
print('interactions_left:',len(F1))
#Filter 2 : only keeping interactions where reference/canonical isoforms had interaction positive
# a two column table which contains refernce isoform which interacts with h.oreofome proteins(proteins in interactor table)
ref_int_pos=F1[(F1['Category']=='reference') & (F1['Interaction_Found']=='positive')][['Gene_Symbol','Interactor_Symbol']].drop_duplicates()
# a file that contains only alternative isoforms
alt_isoforms=F1[F1['Category']=='alternative'].copy()
# adds new column in alt_isoforms which tells wether gene-interactor pair is in ref_int_positive
alt_isoforms['ref_int_was_pos']=alt_isoforms.apply(
    lambda row :((ref_int_pos['Gene_Symbol']==row['Gene_Symbol'])& (ref_int_pos['Interactor_Symbol']==row['Interactor_Symbol'])).any(),
    axis=1)
#only the rows where ref_int_was_pos is True are kept
F2=alt_isoforms[alt_isoforms['ref_int_was_pos']==True]
print('Only_ref_pos:',len(F2))
#Filter 3 : removing genes which are not present in protein_predictions
genes_in_protein_predictions=set(protein_predictions['prot_id'].str.replace(r'_\d+$','',regex=True))
F3=F2[F2['Gene_Symbol'].isin(genes_in_protein_predictions)].copy()
print('genes not in protein predictions :',len(F3))
missing_genes=F2[~F2['Gene_Symbol'].isin(genes_in_protein_predictions)]['Gene_Symbol'].unique()
print('genes removed:',missing_genes)
#Saving all filtered files
F3.to_csv('Filtered_interactor_table.csv',index=False)
protein_predictions.to_csv('ensemble_1423.csv',index=False)
print('Complete')