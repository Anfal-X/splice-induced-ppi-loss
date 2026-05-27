import requests # requests is a library for making HTTP calls (like visiting a website, but from Python).
import time # used to make the script "wait" between checks — so the server is not spammed.
import json # handles JSON format, which is how the API sends back data.
from Bio import SeqIO

fasta_file = "canonical_sequences_filtered.fasta"
sequences = {} #gene_name → sequence string
for record in SeqIO.parse(fasta_file, "fasta"):
    sequences[record.id] = str(record.seq)
    print(f"Loaded {len(sequences)} sequences to submit.")

# Set up the API settings 
BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
EMAIL = "anfalnabeel3@gmail.com"
POLL_INTERVAL = 30 # 30 seconds to wait between checking if the job is done.

# Function to submit ONE sequence 
def submit_sequence(gene_name, sequence):
    url = f"{BASE_URL}/run"
    data = {
        "email": EMAIL,
        "title": gene_name,
        "sequence": sequence,
        "database": "pfam",
        "stype": "p",# "p" means protein sequence
    }

    response = requests.post(url, data=data) # requests.post() sends a POST request to the URL with our data.
    if response.status_code == 200:# status_code 200 means "OK" — the request was accepted successfully.
        job_id = response.text.strip() # response.text is the raw text the server sent back.
        print(f"  Submitted {gene_name} → job ID: {job_id}")
        return job_id# Return the job ID so it can be used later to check the result.

    else:
        print(f"  ERROR submitting {gene_name}: {response.status_code} {response.text}")
        return None

# Function to check if a job is finished 
def check_status(job_id):
    url = f"{BASE_URL}/status/{job_id}"
    response = requests.get(url)# requests.get() fetches a page — like typing a URL in your browser.
    return response.text.strip() # Returns the status text, e.g. "RUNNING" or "FINISHED".

#  Function to download results for a finished job 
def get_results(job_id):
    url = f"{BASE_URL}/result/{job_id}/tsv"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text # Return the raw TSV text. We'll save it to a file in the next step.
    else:
        print(f"  ERROR fetching results for {job_id}: {response.status_code}")
        return None

# Main loop — submit all sequences, wait, collect results 

all_results = []
job_map = {} #mapping job_id → gene_name.

#  Submit all sequences
print("\n--- Submitting all sequences ---")

for gene_name, sequence in sequences.items():
    job_id = submit_sequence(gene_name, sequence)
    if job_id:
        job_map[job_id] = gene_name # If submission succeeded (job_id is not None), store it in our map.
    time.sleep(2)

print(f"\nSubmitted {len(job_map)} jobs. Now waiting for results...")

# Poll until all jobs are done
print("\n--- Checking job statuses ---")

pending_jobs = dict(job_map) # Copy job_map into pending_jobs. 
while pending_jobs:
    finished = []
    for job_id, gene_name in list(pending_jobs.items()):
        status = check_status(job_id)
        print(f"  {gene_name}: {status}")

        if status == "FINISHED":
            result_tsv = get_results(job_id)

            if result_tsv:
                all_results.append((gene_name,result_tsv))

            finished.append(job_id)

        elif status == "ERROR" or status == "FAILED":
            print(f"  !! Job failed for {gene_name} — skipping.")
            finished.append(job_id)
            
    for job_id in finished:
        del pending_jobs[job_id]
       

    if pending_jobs:
        print(f"\n  {len(pending_jobs)} jobs still running. Waiting {POLL_INTERVAL}s...\n")
        time.sleep(POLL_INTERVAL)
        
#  Save all results to one TSV file 

output_file = "interpro_annotations.tsv"

with open(output_file, "w") as f:

    header = "Gene_Symbol\tDB\tDomain_ID\tDomain_name\tStart\tEnd\tScore\tStatus\n"
    f.write(header)

    for gene_name, result_tsv in all_results:
        for line in result_tsv.strip().split("\n"):

            if not line or line.startswith("#"):
                continue

            columns = line.split("\t")
            
            if len(columns) < 10:
                continue
                # Safety check 

            db = columns[3]
            domain_id = columns[4]
            domain_name = columns[5]
            start = columns[6]
            end = columns[7]
            score = columns[8]
            status = columns[9]
            

            if status != "T":
                continue

            f.write(f"{gene_name}\t{db}\t{domain_id}\t{domain_name}\t{start}\t{end}\t{score}\tT\n")

print(f"\nDone! Results saved to: {output_file}")