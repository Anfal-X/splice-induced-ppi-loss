import pandas as pd
File='updated-web-dnet-corrected.csv'
pipenn_results=pd.read_csv(File)
filtered_table=pd.read_csv('Filtered_interactor_table.csv')
alignment_results=pd.read_csv('spliced_positions_reference1.csv')
print('P_R:',len(pipenn_results),
      'F_T:',len(filtered_table),
      'A_R:',len(alignment_results))
# cleaning of pipenn results
pipenn_results['sequence']=pipenn_results['prot_seq'].str.replace(',','',regex=False)# removes commas from the sequences
print(pipenn_results['y_preds'][0][:80])#print first 80 characters of the predictions of the first protein sequence
pipenn_results['y_preds_list']=pipenn_results['y_preds'].str.split(',')# one string is converted into list of individual strings
# every individual string is being converted into actual number
pipenn_results['y_preds_list']=pipenn_results['y_preds_list'].apply(lambda x:# for each row call it x perform the following on it
                                                                    [float(i)for i in x if i.strip()!=''])# ptotects against empty entries
print(pipenn_results['y_preds_list'][0][:5])# prints first 5 scores in the list 
Final_ensemble_table=[]
# looking into filtered table
for row_number,data_in_row in filtered_table.iterrows():# takes out the following for every row(isoform)
    Isoform_id=data_in_row['Isoform_ID']
    Gene_Symbol=data_in_row['Gene_Symbol']
    Interactor_Symbol=data_in_row['Interactor_Symbol']
    Interaction_outcome=data_in_row['Interaction_Found']
    # looking into alignment results
    spliced_row=alignment_results[alignment_results['Isoform_ID']==Isoform_id]
    if len(spliced_row)==0:# no matches found in alignment results
        continue
    spliced_positions=spliced_row['Spliced_positions'].values[0]# taking out of spliced position values from the row
    spliced_positions=eval(spliced_positions)# converting the string into Python list of numbers
    # looking into pipenn results
    pipenn_results_row=pipenn_results[pipenn_results['prot_id']==Isoform_id]
    if len(pipenn_results_row)==0:
        continue
    predictions=pipenn_results_row['y_preds_list'].values[0]# takes out list of predictions
    # goes over every positions and takes out its scores also performs security check so the length of predictions is not exceeded
    spliced_prediction_scores=[predictions[i] for i in spliced_positions if i<len(predictions)]
    # dictionary that saves everything in until the loop ends
    Final_ensemble_table.append({
        'Isoform_ID':Isoform_id,
        'Gene_Symbol':Gene_Symbol,
        'Interactor_Symbol':Interactor_Symbol,
        'Interaction_Outcome':Interaction_outcome,
        'Spliced_Positions': spliced_positions,
        'Number_of_spliced_positions':len(spliced_positions),
        'PiPENN_Prediction_Scores':spliced_prediction_scores,
        'Mean_Scores':sum(spliced_prediction_scores)/len(spliced_prediction_scores) if len(spliced_prediction_scores)>0 else None,
        'Full_PiPENN_Scores':predictions
    })

Final_ensemble_table_df=pd.DataFrame(Final_ensemble_table)
print('Rows:',len(Final_ensemble_table_df))
print(Final_ensemble_table_df.head())
File_name='preds-ensemble-all_final_table_reference1.csv'
Final_ensemble_table_df.to_csv(File_name, index=False)# row numbers are not saved as columns
print('Saved as:',File_name)