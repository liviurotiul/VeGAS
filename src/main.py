import argparse
import subprocess
import os
import sys
import shutil
import warnings
import re

from tqdm import tqdm
from glob import glob

warnings.filterwarnings("ignore")

# ================= Argument Parsing =================
def parse_arguments():
    parser = argparse.ArgumentParser(description="Run the PeGAS pipeline.")
    parser.add_argument("-d", "--data", dest="data", help="Directory containing all the fastq.gz files", required=True)
    parser.add_argument("-o", "--output", dest="output", help="Directory where output files will be saved", required=True)
    parser.add_argument("-c", "--cores", dest="cores", help="The number of cores to use", default=4, type=int)
    parser.add_argument("-r", "--reference", dest="reference", help="Path to folder of reference genomes", required=True)
    parser.add_argument("-t", "--host", dest="host", help="Path to host genomes", required=True)
    parser.add_argument("-w", "--overwrite", dest="overwrite", help="Overwrite the output directory if it exists", action="store_true")
    parser.add_argument("-cc", "--context_cores", dest="context_cores", help="The number of cores to use for context", default=1, type=int)
    return parser.parse_args()

# ================= Utility Functions =================
def list_fastq_files(path):
    """Returns a list of all .fastq.gz files in the specified path relative to base_folder."""
    full_path = os.path.join(path)
    return [f for f in glob(os.path.join(full_path, "*.fastq.gz"))]

def get_core_sample_name(filename):
    """Extracts the core sample name by removing _R1 or _R2 and other suffixes."""
    return os.path.basename(filename).replace("_R1", "").replace("_R2", "").replace(".fastq.gz", "")

def build_fastq_pairs(fastq_files):
    """Pairs R1 and R2 files based on sample names."""
    pairs = {}
    for file in fastq_files:
        sample = get_core_sample_name(file)
        if sample:
            if sample not in pairs:
                pairs[sample] = {}
            if "_R1" in file or "R1" in file:
                pairs[sample]["R1"] = file
            elif "_R2" in file or "R2" in file:
                pairs[sample]["R2"] = file
    return {s: p for s, p in pairs.items() if "R1" in p and "R2" in p}

def copy_files(base_folder, source_files, destination="raw_data"):
    """Copies files from source directory to destination within base_folder without renaming them."""
    dest_path = os.path.join(base_folder, destination)
    os.makedirs(dest_path, exist_ok=True)
    for src in tqdm(source_files, desc="Copying files"):
        dest = os.path.join(dest_path, os.path.basename(src))
        
        # Move _R1 and _R2 to the end of the filename
        if "R1" in dest:
            dest = dest.replace("_R1", "").replace(".fastq.gz", "_R1.fastq.gz")
        elif "R2" in dest:
            dest = dest.replace("_R2", "").replace(".fastq.gz", "_R2.fastq.gz")
        
        if not os.path.exists(dest) or os.path.getsize(src) != os.path.getsize(dest):
            shutil.copy2(src, dest)
            tqdm.write(f"[pegas] Copied '{src}' to '{dest}'.")

def remove_extra_files(base_folder, destination, valid_files):
    """Removes unwanted files from the destination directory and clears related data for affected samples."""
    dest_path = os.path.join(base_folder, destination)
    valid_basenames = {os.path.basename(f) for f in valid_files}
    valid_samples = list(set([get_core_sample_name(f) for f in valid_files]))

    # # Remove extra files in the destination directory
    # for file in glob(os.path.join(dest_path, "*")):
    #     if os.path.basename(file) not in valid_basenames:
    #         os.remove(file)
    #         tqdm.write(f"[pegas] Removed '{file}'.")

    # # Remove extra sample directories in the results folder
    # results_path = os.path.join(base_folder, "results")
    # if os.path.exists(results_path):
    #     for folder in os.listdir(results_path):
    #         if folder not in valid_samples:
    #             shutil.rmtree(os.path.join(results_path, folder))
    #             tqdm.write(f"[pegas] Removed 'results/{folder}' directory.")

    # # Remove outdated files in fastqc directory
    # fastqc_path = os.path.join(base_folder, "fastqc")
    # if os.path.exists(fastqc_path):
    #     valid_file_names = [os.path.basename(file).replace(".fastq.gz", "") for file in valid_files]
    #     for file in os.listdir(fastqc_path):
    #         if os.path.basename(file).replace("_fastqc.html", "").replace("_fastqc.zip", "") not in valid_file_names:
    #             os.remove(os.path.join(fastqc_path, file))
    #             tqdm.write(f"[pegas] Removed '{file}'.")

def main():

    install_path = os.path.dirname(os.path.realpath(__file__))

    args = parse_arguments()

    data_dir = args.data
    output_dir = args.output
    cores = args.cores
    overwrite = args.overwrite
    reference = args.reference
    host = args.host
    ccores = args.context_cores
    
    # Check if the data directory exists
    if not os.path.exists(data_dir):
        tqdm.write(f"[pegas]Data directory '{data_dir}' does not exist.")
        sys.exit(1)
    if os.listdir(data_dir) == []:
        tqdm.write(f"[pegas]Data directory '{data_dir}' is empty.")
        sys.exit(1)
    if not os.path.exists(reference):
        tqdm.write(f"[pegas]Reference directory '{reference}' does not exist.")
        sys.exit(1)
    if os.listdir(reference) == []:
        tqdm.write(f"[pegas]Reference directory '{reference}' is empty.")
        sys.exit(1)
    if not os.path.exists(host):
        tqdm.write(f"[pegas]Host directory '{host}' does not exist.")
        sys.exit(1)
    if os.listdir(host) == []:
        tqdm.write(f"[pegas]Host directory '{host}' is empty.")
        sys.exit(1)
    
    # Build config params as a list (order preserved) to pass to snakemake
    config_params = [
        f"install_path={install_path}",
        f"output_dir={output_dir}",
        f"host_genome={host}",
        f"reference_genome={reference}",
        f"ccores={ccores}",
    ]
    
    # List all FASTQ files in the raw_data_path and raw_data directories
    raw_data_files = list_fastq_files(data_dir)

    # Check if the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        tqdm.write(f"[pegas]Created output directory '{output_dir}'")
    else:
        tqdm.write(f"[pegas]Using existing output directory '{output_dir}'")

    # Copy new or modified files from raw_data_path to raw_data
    copy_files(output_dir, raw_data_files, "raw_data")

    # Remove files in output that do not exist in raw_data_path
    remove_extra_files(output_dir, "raw_data", raw_data_files)

    # Build the Snakemake command
    unlock_command = [
        "snakemake",
        "--snakefile", os.path.join(install_path, "Snakefile"),
        "--directory", output_dir,
        "--cores", str(cores),
        "--rerun-incomplete",
        "--use-conda",
        "--unlock"
    ]
    
    command = [
        "snakemake",
        "--snakefile", os.path.join(install_path, "Snakefile"),
        "--directory", output_dir,
        "--cores", str(cores),
        "--rerun-incomplete",
        "--use-conda"
    ]

    command.append("--config")
    unlock_command.append("--config")
    command.extend(config_params)
    unlock_command.extend(config_params)

    # Run the pipeline
    print("Running command:" + " ".join(command))
    # print("Unlock command:" + " ".join(unlock_command))
    # Try unlocking first (non-fatal if it fails)
    try:
        subprocess.run(unlock_command, check=False)
    except Exception:
        tqdm.write("[pegas] Warning: unlock command failed (continuing)")

    result = subprocess.run(command)
    if result.returncode != 0:
        tqdm.write("Error: Pipeline failed.")
        sys.exit(result.returncode)
    else:
        tqdm.write("Pipeline completed successfully.")
        igv_report_path = os.path.abspath(os.path.join(output_dir, "igv_report.html"))
        tqdm.write(f"\n========================================")
        tqdm.write(f"Raportul IGV HTML este gata!")
        tqdm.write(f"Îl poți accesa din browser la adresa:")
        tqdm.write(f"file://{igv_report_path}")
        tqdm.write(f"========================================\n")

if __name__ == "__main__":
    main()
