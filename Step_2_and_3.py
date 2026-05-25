
import pandas as pd
import ast
import os
import glob
DATA_FOLDER = "/Users/user/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Internship/Code"
DOMAINS_FILE = "interpro_pfam_only.tsv"
OUTPUT_FOLDER = "domain_analysis_results"
PIPENN_THRESHOLD = 0.5
MODEL_NAMES = ["ensemble", "ann", "rnn", "cnn-rnn", "dnet", "web-dnet", "rnet", "unet"]

# Creates output folder structure 
os.makedirs(f"{OUTPUT_FOLDER}/group1_reference", exist_ok=True)#creates folder and subfolder per group
os.makedirs(f"{OUTPUT_FOLDER}/group2_longest", exist_ok=True)
os.makedirs(f"{OUTPUT_FOLDER}/MO-GI", exist_ok=True)# exist_ok=True means "don't crash if the folder already exists".
print("Output folders created.")

# Detects which group a file belongs to 

def detect_group(filename):
    name = os.path.basename(filename).lower()
    # .lower() makes matching case-insensitive.

    # Detect model name first.
    detected_model = None
    for model in sorted(MODEL_NAMES, key=len, reverse=True):
        if model.replace("-", "") in name.replace("-", ""):
            detected_model = model
            break

    
    if "triples" in name:
        return "MO-GI", detected_model

    elif "longest" in name:
        return "group2_longest", detected_model

    elif "reference" in name:
        return "group1_reference", detected_model

    else:
        return None, None
    # We moved "triple" to the top because the previous version
    # was accidentally matching something else in the filename first.
# Start of Step 3:Expand + merge domains
# Helper functions 

def parse_positions(value):
    # Converts "[32, 33, 34]" string → real Python list [32, 33, 34]
    try:
        result = ast.literal_eval(value)
        if isinstance(result, list):
            return [int(x) for x in result]
        else:
            return [int(result)]
    except:
        return []

def parse_scores(value):
    # Converts "[0.8, 0.9, 0.7]" string → real Python list [0.8, 0.9, 0.7]
    try:
        result = ast.literal_eval(value)
        if isinstance(result, list):
            return [float(x) for x in result]
        else:
            return [float(result)]
    except:
        return []

def expand_isoforms(df):
    # Takes final table and expands it to one row per spliced residue.

    df["Spliced_Positions_parsed"] = df["Spliced_Positions"].apply(parse_positions)
    df["Scores_parsed"] = df["PiPENN_Prediction_Scores"].apply(parse_scores)

    rows = []
    for _, row in df.iterrows():
        positions = row["Spliced_Positions_parsed"]
        scores    = row["Scores_parsed"]

        for i, pos in enumerate(positions):
            score = scores[i] if i < len(scores) else None
            rows.append({
                "Isoform_ID":          row["Isoform_ID"],
                "Gene_Symbol":         row["Gene_Symbol"],
                "Interactor_Symbol":   row["Interactor_Symbol"],
                "Interaction_Outcome": row["Interaction_Outcome"],
                "Residue_Position":    pos,
                "PiPENN_Score":        score,
                "Is_Interface":        score >= PIPENN_THRESHOLD if score is not None else False,
                "Mean_Score":          row["Mean_Scores"],
                "Number_Spliced":      row["Number_of_spliced_positions"],
            })

    return pd.DataFrame(rows)

def expand_domains(domains_df):
    # Expands domain start-end ranges to one row per residue.

    domain_rows = []
    for _, row in domains_df.iterrows():
        start = int(row["Start"])
        end   = int(row["End"])
        for pos in range(start, end + 1):
            domain_rows.append({
                "Gene_Symbol":      row["Gene_Symbol"],
                "Residue_Position": pos,
                "Domain_ID":        row["Domain_ID"],
                "Domain_name":      row["Domain_name"],
                "Domain_Start":     start,
                "Domain_End":       end,
                "Domain_Length":    end - start + 1,
            })
    return pd.DataFrame(domain_rows)

# Start of Step 4: Calculation of prop_A/B/C
def calculate_proportions(merged_df):
    # calculates prop_A, prop_B, prop_C per isoform × domain.
    # merged_df is the joined table of residues + domain annotations.

    results = []

    # Group by every unique combination of isoform + domain.
    grouped = merged_df[merged_df["In_Domain"]].groupby(
        ["Isoform_ID", "Gene_Symbol", "Interactor_Symbol",
         "Interaction_Outcome", "Domain_ID", "Domain_name",
         "Domain_Start", "Domain_End", "Domain_Length"]
    )

    for keys, group in grouped:
        # keys = the values of the groupby columns for this group (as a tuple)
        # group = a mini-table of all rows in this group

        isoform_id, gene, interactor, outcome, dom_id, dom_name, dom_start, dom_end, dom_len = keys
        dom_start = int(dom_start)
        dom_end   = int(dom_end)
        dom_len   = int(dom_len)

        all_spliced_in_isoform = merged_df[
            merged_df["Isoform_ID"] == isoform_id
        ]["Residue_Position"].unique()
        # Filter merged_df to only rows for this isoform.
        # Get all unique residue positions that are spliced in this isoform.

        # Residues in this domain that are also spliced in this isoform.
        domain_positions = set(range(dom_start, dom_end + 1))

        spliced_positions = set(all_spliced_in_isoform)
        # All spliced residue positions for this isoform, as a set.

        interface_positions = set(
            merged_df[
                (merged_df["Isoform_ID"] == isoform_id) &
                (merged_df["Is_Interface"] == True)
            ]["Residue_Position"].unique()
        )

        #  The 3 key proportions 

        prop_A = len(domain_positions & interface_positions) / dom_len
        # domain_positions & interface_positions = residues that are BOTH:
        #   - inside this domain
        #   - predicted as interface
        # Divided by domain length = proportion of the domain that is interface.

        prop_B = len(domain_positions & spliced_positions) / dom_len
        # domain_positions & spliced_positions = residues that are BOTH:
        #   - inside this domain
        #   - spliced out in this isoform
        
        prop_C = len(domain_positions & spliced_positions & interface_positions) / dom_len
        # MO-GI intersection: residues that are ALL THREE:
        #   - inside this domain
        #   - spliced out
        #   - predicted as interface
        
        if prop_B == 0 and prop_C == 0:
            continue
            # Skip this combination if no splicing occurs in this domain.
            

        results.append({
            "Isoform_ID":          isoform_id,
            "Gene_Symbol":         gene,
            "Interactor_Symbol":   interactor,
            "Interaction_Outcome": outcome,
            "Domain_ID":           dom_id,
            "Domain_name":         dom_name,
            "Domain_Start":        dom_start,
            "Domain_End":          dom_end,
            "Domain_Length":       dom_len,
            "prop_A":              round(prop_A, 4),
            "prop_B":              round(prop_B, 4),
            "prop_C":              round(prop_C, 4),
            
        })

    return pd.DataFrame(results)

#  Load domains ONCE (shared across all 24 datasets) 

print("Loading Pfam domain annotations...")
domains_df = pd.read_csv(DOMAINS_FILE, sep="\t")
domain_residues = expand_domains(domains_df)
print(f"Domain residues expanded: {len(domain_residues)} rows")

# Find all 24 CSV files and process them 

all_files = glob.glob(f"{DATA_FOLDER}/*.csv")
# Returns a list of full file paths.

print(f"\nFound {len(all_files)} CSV files in folder.")

summary = []

for filepath in sorted(all_files):#  processes files in alphabetical order 

    filename = os.path.basename(filepath)
    group, model = detect_group(filepath)
    
    if group is None or model is None:
        print(f"  SKIPPED (could not classify): {filename}")
        continue
        
    print(f"\nProcessing: {filename}")
    print(f"  Group: {group} | Model: {model}")

    # Load the dataset 
    df = pd.read_csv(filepath)

    #  Expand isoforms to residue level and merge with domains 
    residue_df = expand_isoforms(df)
    # Expand each isoform row into one row per spliced residue.

    merged = pd.merge(
        residue_df,
        domain_residues,
        on=["Gene_Symbol", "Residue_Position"],
        how="left"
    )
    

    merged["In_Domain"] = merged["Domain_ID"].notna()
    # True if the residue falls inside a Pfam domain, False otherwise.

    #  Calculate prop_A, prop_B, prop_C 
    proportions_df = calculate_proportions(merged)

    # Save results 
    out_path = f"{OUTPUT_FOLDER}/{group}/{model}_proportions.csv"
    proportions_df.to_csv(out_path, index=False)
    

    print(f"  Isoform-residues: {len(residue_df)} | In-domain: {merged['In_Domain'].sum()} | Prop rows: {len(proportions_df)}")
    print(f"  Saved → {out_path}")

    summary.append({
        "File":       filename,
        "Group":      group,
        "Model":      model,
        "Isoforms":   len(df),
        "Residues":   len(residue_df),
        "Prop_rows":  len(proportions_df),
    })

# Print final summary 
print("\n" + "="*60)
print("PROCESSING COMPLETE — SUMMARY")
print("="*60)
summary_df = pd.DataFrame(summary)
print(summary_df.to_string(index=False))


summary_df.to_csv(f"{OUTPUT_FOLDER}/processing_summary.csv", index=False)
print(f"\nSummary saved to: {OUTPUT_FOLDER}/processing_summary.csv")