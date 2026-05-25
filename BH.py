
import pandas as pd
from scipy import stats
import numpy as np

File = pd.read_csv('all_pvalues.csv')

print("File loaded successfully!")
print(f"Number of rows: {len(File)}")

wilcoxon_raw = File['Wilcoxon_p'].values # .values converts it from a pandas column into a plain numpy array
jaccard_raw = File['Jaccard_p'].values
chi2_raw = File['Chi2_p'].values

wilcoxon_bh = stats.false_discovery_control(wilcoxon_raw, method='bh')
jaccard_bh = stats.false_discovery_control(jaccard_raw, method='bh')
chi2_bh = stats.false_discovery_control(chi2_raw, method='bh')

def significance_marker(p):
  
    if p < 0.001:
        return '***'   # very highly significant
    elif p < 0.01:
        return '**'    # highly significant
    elif p < 0.05:
        return '*'     # significant
    else:
        return 'ns'    # not significant (ns = not significant)


results = pd.DataFrame()
results['Model'] = File['File']
results['Wilcoxon_raw_p'] = wilcoxon_raw
results['Wilcoxon_corrected_p'] = wilcoxon_bh
results['Wilcoxon_significance'] = [significance_marker(p) for p in wilcoxon_bh]

results['Jaccard_raw_p'] = jaccard_raw
results['Jaccard_corrected_p'] = jaccard_bh
results['Jaccard_significance'] = [significance_marker(p) for p in jaccard_bh]

results['Chi2_raw_p'] = chi2_raw
results['Chi2_corrected_p'] = chi2_bh
results['Chi2_significance'] = [significance_marker(p) for p in chi2_bh]

print("\n--- Significance summary (corrected p < 0.05) ---")
for test, col in [('Wilcoxon', 'Wilcoxon_corrected_p'),
                  ('Jaccard',  'Jaccard_corrected_p'),
                  ('Chi2',     'Chi2_corrected_p')]:
    # this loop runs 3 times, once for each test
    # each time, 'test' is the test name and 'col' is the column name
    n_sig = (results[col] < 0.05).sum()
    # results[col] < 0.05 creates a True/False list (True where p < 0.05)
    # .sum() counts the True values (True = 1, False = 0)
    print(f"{test}: {n_sig} / {len(results)} significant")


summary = pd.DataFrame({
    'Test':                            ['Wilcoxon', 'Jaccard', 'Chi-square'],
    'Significant raw (p < 0.05)':      [(wilcoxon_raw < 0.05).sum(),
                                        (jaccard_raw  < 0.05).sum(),
                                        (chi2_raw     < 0.05).sum()],
    'Significant BH (p < 0.05)':       [(wilcoxon_bh < 0.05).sum(),
                                        (jaccard_bh  < 0.05).sum(),
                                        (chi2_bh     < 0.05).sum()],
})

output_file = 'BH_corrected_results.xlsx'
# this is the name of the output file that will be created

with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
    

    results.to_excel(writer, sheet_name='BH Corrected Results', index=False)
   

    summary.to_excel(writer, sheet_name='Method Comparison', index=False)
    
print(f"\nDone! Results saved to: {output_file}")

