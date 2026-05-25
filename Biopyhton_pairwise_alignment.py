from Bio.Align import PairwiseAligner
import pandas as pd
#loading of files
pipenn_results=pd.read_csv('updated-web-dnet-corrected.csv')
filtered_table=pd.read_csv('Filtered_interactor_table.csv')
print('number of proteins:',len(pipenn_results))
print('number of interactions:',len(filtered_table))
#omitting commas from the protein sequences
pipenn_results['sequence']=pipenn_results['prot_seq'].str.replace(',','',regex=False)#treat comma as plain character
print('example:',pipenn_results['prot_id'][0])
print('example sequence:',pipenn_results['sequence'][0][:20]) # first 20 sequences 
# takes out unique alternative isoform sequences
isoforms_aligned=filtered_table['Isoform_ID'].unique()
print('alt_isoforms:',len(isoforms_aligned))
# setting up the paramaters of aligner
aligner=PairwiseAligner()
aligner.mode='global'
aligner.match_score=1
aligner.mismatch_score=-1
aligner.open_gap_score=-10
aligner.extend_gap_score=-0.5
# store results 
alignment_results=[]
# Main loop which loops through every alternative isoform
for isoform_id in isoforms_aligned:
    gene_name=isoform_id.rsplit('_',1)[0]#removes number from the end of gene symbols
    gene_rows=pipenn_results[pipenn_results['prot_id'].str.startswith(gene_name+'_')]
    if len(gene_rows)==0:
        print('No sequences found for gene:',gene_name)
        continue
    longest_row=gene_rows.loc[gene_rows['sequence'].str.len().idxmax()]#finds the longest sequence for the gene
    longest_id=longest_row['prot_id']
    if longest_id==gene_name+'_1':
        reference_id=gene_name+'_1'
        reference_type='_1'
    else:
        reference_id=longest_id
        reference_type='longest'
    ref_row=pipenn_results[pipenn_results['prot_id']==reference_id]# results in reference sequence 
    if len(ref_row)==0:
        print('No reference id:',isoform_id)
        continue
    ref_seq= ref_row['sequence'].values[0]
    alt_row=pipenn_results[pipenn_results['prot_id']==isoform_id]# gives alternative isoform
    if len(alt_row)==0:
        print('No alternative isoform:',isoform_id)
        continue
    alt_seq=alt_row['sequence'].values[0]
    # run the alignment
    run_alignment=aligner.align(ref_seq,alt_seq)
    best_alignment=run_alignment[0]
    aligned_ref=best_alignment.aligned[0]
    aligned_alt=best_alignment.aligned[1]
    # set of the positions of reference sequence that are part of alternative sequence
    ref_part_alt=set()
    for ref_block,alt_block in zip(aligned_ref,aligned_alt):
        ref_start=ref_block[0]
        ref_end=ref_block[1]
        for pos in range (ref_start,ref_end):
            ref_part_alt.add(pos)
    # taking out spliced out regions
    all_ref_positions=set(range(len(ref_seq)))
    spliced_out_positions=sorted(all_ref_positions-ref_part_alt)
    if len(spliced_out_positions)==0:# removing of true gain of sequence isoforms
        print('Excluded (gain-of-sequence):',isoform_id)
        continue
    # store the results
    alignment_results.append({
        'Isoform_ID': isoform_id,
        'Gene': gene_name,
        'Reference_ID': reference_id,
        'Reference_type': reference_type,
        'Ref_length': len(ref_seq),
        'Alt_length': len(alt_seq),
        'Spliced_positions': spliced_out_positions,
        'Num_spliced': len(spliced_out_positions)
    })
    print('Done:',isoform_id,'spliced out positions:',len(spliced_out_positions))
# Saving results
alignment_results_df = pd.DataFrame(alignment_results)

reference_1_df = alignment_results_df[alignment_results_df['Reference_type'] == '_1']
reference_longest_df = alignment_results_df[alignment_results_df['Reference_type'] == 'longest']

reference_1_df.to_csv('spliced_positions_reference1.csv', index=False)
reference_longest_df.to_csv('spliced_positions_longest.csv', index=False)

print('Reference _1 isoforms saved:', len(reference_1_df))
print('Longest reference isoforms saved:', len(reference_longest_df))
print('Total saved:', len(alignment_results_df))

