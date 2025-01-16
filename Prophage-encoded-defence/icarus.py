#icarus adjusted for pf candidates, removing the requirement for DUF2130

import os
import pandas as pd
import requests
import subprocess
import time

API_KEY = "ad86a968909287348cb03552eff8dbde7b08"  # Replace with your NCBI API key or set to None if not using
DELAY = 5  # Longer delay (seconds) between requests to avoid bans

def fetch_sequences(wp_ids, output_file):
    """Fetch FASTA sequences for a list of WP IDs and save to a file."""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    batch_size = 10  # Split requests into smaller batches
    with open(output_file, "w") as f:
        for i in range(0, len(wp_ids), batch_size):
            batch = wp_ids[i:i + batch_size]
            params = {
                "db": "protein",
                "rettype": "fasta",
                "id": ",".join(batch),
            }
            if API_KEY:
                params["api_key"] = API_KEY
            response = requests.get(base_url, params=params)
            if response.status_code == 200:
                f.write(response.text)
                time.sleep(DELAY)  # Longer delay to avoid bans
            elif response.status_code == 429:
                print("Rate limit exceeded. Retrying after a longer delay...")
                time.sleep(DELAY * 2)  # Double the delay for retry
                i -= batch_size  # Retry the same batch
            else:
                raise Exception(f"Failed to fetch sequences. HTTP status: {response.status_code}")

def find_dominant_protein(fasta_file, keyword):
    """Identify a protein in the FASTA file with the given keyword in its header."""
    dominant_protein = None
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">") and keyword in line:
                dominant_protein = line.split()[0][1:]  # Extract protein ID from header
                break
    if not dominant_protein:
        raise Exception(f"No protein with keyword '{keyword}' found in {fasta_file}.")
    return dominant_protein

def run_cblaster(input_fasta, output_dir, comment, num_genes):
    """Run cblaster with the given input FASTA and save results to the output directory."""
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{comment}_cblaster_results.json")
    binary_file = os.path.join(output_dir, f"{comment}_cblaster_results.binary")
    cblaster_cmd = [
        "cblaster",
        "search",
        "-qf", input_fasta,
        "-o", output_file,
        "-b", binary_file,
        "-mi", "20",
        "-mc", "20",
        "-mh", str(num_genes),  # Set -mh to the number of genes in the operon
        "-g", "70",
        # "-r", dominant_protein,  # Commented out -r flag
        "-hs", "500"
    ]
    subprocess.run(cblaster_cmd, check=True)

def skip_if_processed(output_dir):
    """Check if output directory exists and contains files."""
    if os.path.exists(output_dir) and any(os.scandir(output_dir)):
        print(f"Skipping {output_dir} as it already contains files.")
        return True
    return False

def main():
    # Load the table into a Pandas DataFrame
    table_file = "pf-operons.txt"  # Replace with your table file path
    df = pd.read_csv(table_file, sep="\t")

    # Group by nuccore_id and Comment
    grouped = df.groupby(["nuccore_id", "Comment"])

    for (nuccore_id, comment), group in grouped:
        if pd.isna(comment):
            continue  # Skip rows with no comment

        # Define file paths
        output_dir = comment.replace(" ", "_")
        fasta_file = f"{nuccore_id}_{comment}_sequences.fasta"

        # Skip if already processed
        if skip_if_processed(output_dir):
            continue

        # Get WP IDs from the group
        wp_ids = group["Protein_ID"].tolist()

        # Ensure clean slate for each operon
        if os.path.exists(fasta_file):
            os.remove(fasta_file)

        # Fetch FASTA sequences
        print(f"Fetching sequences for WP IDs: {wp_ids}")
        try:
            fetch_sequences(wp_ids, fasta_file)
        except Exception as e:
            print(f"Error processing operon {comment}: {e}")
            continue  # Skip to the next operon

        # Calculate the number of genes in the operon
        num_genes = len(wp_ids)

        # Commented out dominant protein search for "DUF2130"
        # print(f"Searching for dominant protein with keyword 'DUF2130' in {fasta_file}")
        # try:
        #     dominant_protein = find_dominant_protein(fasta_file, "DUF2130")
        #     print(f"Dominant protein found: {dominant_protein}")
        # except Exception as e:
        #     print(f"Error finding dominant protein: {e}")
        #     continue  # Skip to the next operon

        # Run cblaster
        print(f"Running cblaster for comment: {comment}")
        try:
            run_cblaster(fasta_file, output_dir, comment, num_genes)
        except subprocess.CalledProcessError as e:
            print(f"Error running cblaster: {e}")
            continue  # Skip to the next operon

        # Clean up FASTA file after processing
        if os.path.exists(fasta_file):
            os.remove(fasta_file)

        print(f"Completed processing for operon: {comment}\n")

if __name__ == "__main__":
    main()
