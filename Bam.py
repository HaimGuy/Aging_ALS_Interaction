import os
import subprocess
import pandas as pd
import tkinter as tk
from tkinter import filedialog, messagebox
from rpy2.robjects import pandas2ri
from rpy2.robjects import r
from rpy2.robjects.packages import importr

# Activate pandas2ri for R-Python conversion
pandas2ri.activate()

# Import DESeq2 from R
base = importr("base")
utils = importr("utils")
DESeq2 = importr("DESeq2")

# Utility Functions
def generate_counts_matrix(bam_folder, annotation_file, output_folder):
    """
    Generate counts matrix from BAM files using featureCounts.
    """
    os.makedirs(output_folder, exist_ok=True)
    bam_files = [os.path.join(bam_folder, f) for f in os.listdir(bam_folder) if f.endswith(".bam")]
    
    counts_file = os.path.join(output_folder, "counts.txt")
    command = [
        "featureCounts",
        "-a", annotation_file,
        "-o", counts_file,
        "-T", "4",  # Number of threads
        "-s", "0",  # Unstranded data
    ] + bam_files

    print(f"Running featureCounts...")
    subprocess.run(command, check=True)
    print(f"Counts matrix saved to {counts_file}")
    return counts_file

def filter_low_counts(counts_df, threshold):
    """
    Filters out genes with total counts below the specified threshold.
    """
    total_counts = counts_df.sum(axis=1)  # Sum across all samples for each gene
    filtered_counts_df = counts_df[total_counts >= threshold]
    print(f"Filtered genes: {len(counts_df) - len(filtered_counts_df)} genes removed (below threshold of {threshold}).")
    return filtered_counts_df

def perform_differential_expression(counts_file, metadata_file, group_column, output_file, read_count_threshold):
    """
    Perform differential expression analysis using DESeq2 and save results to Excel.
    """
    print("Loading counts and metadata...")
    # Load counts matrix
    counts_df = pd.read_csv(counts_file, sep="\t", comment="#", index_col=0)
    counts_df = counts_df.iloc[:, 5:]  # Remove non-count columns
    
    # Apply read count threshold filter
    counts_df = filter_low_counts(counts_df, read_count_threshold)

    # Load metadata
    metadata_df = pd.read_csv(metadata_file)

    # Prepare R data structures
    r_counts = pandas2ri.py2rpy(counts_df)
    r_metadata = pandas2ri.py2rpy(metadata_df)

    print("Running DESeq2...")
    dds = DESeq2.DESeqDataSetFromMatrix(
        countData=r_counts,
        colData=r_metadata,
        design=r(f"~ {group_column}")
    )
    dds = DESeq2.DESeq(dds)
    res = DESeq2.results(dds)

    print("Processing results...")
    res_df = pandas2ri.rpy2py(res).reset_index()
    res_df.columns = ["Gene", "Log2 Fold Change", "P-value", "Adjusted P-value"]
    
    # Save results to Excel
    res_df.to_excel(output_file, index=False)
    print(f"Differential expression results saved to {output_file}")

# GUI for Input Selection
def start_gui():
    """
    Creates a GUI to select inputs for differential expression analysis.
    """
    root = tk.Tk()
    root.title("Differential Expression Analysis")

    # Variables
    bam_folder = tk.StringVar()
    annotation_file = tk.StringVar()
    metadata_file = tk.StringVar()
    output_file = tk.StringVar()
    group_column = tk.StringVar(value="Group")
    read_count_threshold = tk.IntVar(value=10)  # Default read count threshold

    def select_folder(var):
        folder_path = filedialog.askdirectory(title="Select Folder")
        if folder_path:
            var.set(folder_path)

    def select_file(var):
        file_path = filedialog.askopenfilename(title="Select File")
        if file_path:
            var.set(file_path)

    def save_file(var):
        file_path = filedialog.asksaveasfilename(
            title="Save Output File",
            defaultextension=".xlsx",
            filetypes=[("Excel Files", "*.xlsx"), ("All Files", "*.*")]
        )
        if file_path:
            var.set(file_path)

    def run_analysis():
        try:
            print("Starting analysis...")
            counts_file = generate_counts_matrix(
                bam_folder=bam_folder.get(),
                annotation_file=annotation_file.get(),
                output_folder="counts_output"
            )
            perform_differential_expression(
                counts_file=counts_file,
                metadata_file=metadata_file.get(),
                group_column=group_column.get(),
                output_file=output_file.get(),
                read_count_threshold=read_count_threshold.get()
            )
            messagebox.showinfo("Success", f"Differential expression results saved to {output_file.get()}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    # GUI Layout
    tk.Label(root, text="Select BAM Folder:").pack(pady=5)
    tk.Entry(root, textvariable=bam_folder, width=50).pack(pady=5)
    tk.Button(root, text="Browse", command=lambda: select_folder(bam_folder)).pack(pady=5)

    tk.Label(root, text="Select GTF Annotation File:").pack(pady=5)
    tk.Entry(root, textvariable=annotation_file, width=50).pack(pady=5)
    tk.Button(root, text="Browse", command=lambda: select_file(annotation_file)).pack(pady=5)

    tk.Label(root, text="Select Metadata File (CSV):").pack(pady=5)
    tk.Entry(root, textvariable=metadata_file, width=50).pack(pady=5)
    tk.Button(root, text="Browse", command=lambda: select_file(metadata_file)).pack(pady=5)

    tk.Label(root, text="Group Column Name:").pack(pady=5)
    tk.Entry(root, textvariable=group_column, width=20).pack(pady=5)

    tk.Label(root, text="Read Count Threshold:").pack(pady=5)
    tk.Entry(root, textvariable=read_count_threshold, width=20).pack(pady=5)

    tk.Label(root, text="Save Output Excel File As:").pack(pady=5)
    tk.Entry(root, textvariable=output_file, width=50).pack(pady=5)
    tk.Button(root, text="Save As", command=lambda: save_file(output_file)).pack(pady=5)

    tk.Button(root, text="Run Analysis", command=run_analysis).pack(pady=20)

    root.mainloop()

if __name__ == "__main__":
    start_gui()
