import os
import pandas as pd
import numpy as np
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from scipy.stats import hypergeom
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Utility Functions
def save_excel(dataframes, output_path):
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        for sheet_name, df in dataframes.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    print(f"Results saved to {output_path}")

def plot_venn(aging_genes, als_genes, output_folder):
    overlap = set(aging_genes).intersection(set(als_genes))
    plt.figure()
    venn2([set(aging_genes), set(als_genes)], ("Aging Genes", "ALS Genes"))
    plt.title("Venn Diagram of Gene Intersections")
    venn_path = os.path.join(output_folder, "venn_diagram.png")
    plt.savefig(venn_path)
    plt.show()
    return venn_path

def compare_gene_lists(aging_genes_df, als_genes_df, p_value_threshold, log_fc_threshold):
    """
    Filters and compares Aging and ALS gene lists based on p-value and log fold change thresholds.
    """
    # Filter Aging genes
    aging_filtered = aging_genes_df[
        (aging_genes_df["P-value"] <= p_value_threshold) &
        (abs(aging_genes_df["Log2 Fold Change"]) >= log_fc_threshold)
    ]
    
    # Filter ALS genes
    als_filtered = als_genes_df[
        (als_genes_df["P-value"] <= p_value_threshold) &
        (abs(als_genes_df["Log2 Fold Change"]) >= log_fc_threshold)
    ]
    
    # Find intersecting genes
    intersecting_genes = set(aging_filtered["Gene"]).intersection(set(als_filtered["Gene"]))

    # Extract statistics for intersecting genes
    intersecting_stats = []
    for gene in intersecting_genes:
        # Extract the rows for the current gene from both filtered lists
        aging_row = aging_filtered.loc[aging_filtered["Gene"] == gene].iloc[0]
        als_row = als_filtered.loc[als_filtered["Gene"] == gene].iloc[0]
        intersecting_stats.append({
            "Gene": gene,
            "Aging Log2 Fold Change": aging_row["Log2 Fold Change"],
            "Aging P-value": aging_row["P-value"],
            "ALS Log2 Fold Change": als_row["Log2 Fold Change"],
            "ALS P-value": als_row["P-value"],
        })

    intersecting_stats_df = pd.DataFrame(intersecting_stats)

    # Calculate hypergeometric test
    total_genes = 20000  # Example: total genes in the dataset
    p_overlap = hypergeom.sf(len(intersecting_genes) - 1, total_genes, len(aging_filtered), len(als_filtered))

    return aging_filtered, als_filtered, intersecting_stats_df, p_overlap

# Main Analysis Pipeline
def analyze_data(params, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    dataframes = {}

    # Load Aging and ALS gene data
    aging_genes_df = pd.read_excel(params["aging_file"])
    als_genes_df = pd.read_excel(params["als_file"])

    # Compare the gene lists
    aging_filtered, als_filtered, intersecting_stats, p_overlap = compare_gene_lists(
        aging_genes_df,
        als_genes_df,
        params["p_value_threshold"],
        params["log_fc_threshold"]
    )

    venn_path = plot_venn(aging_filtered["Gene"], als_filtered["Gene"], output_folder)

    # Save results
    dataframes["Aging Filtered Genes"] = aging_filtered
    dataframes["ALS Filtered Genes"] = als_filtered
    dataframes["Intersecting Genes"] = intersecting_stats
    dataframes["Statistics"] = pd.DataFrame([{"Metric": "Overlap P-value", "Value": p_overlap}])

    save_excel(dataframes, os.path.join(output_folder, "analysis_results.xlsx"))

    return venn_path

# GUI
def start_gui():
    root = tk.Tk()
    root.title("ALS Aging Interaction Analysis")

    # Variables for inputs and thresholds
    aging_file = tk.StringVar()
    als_file = tk.StringVar()
    p_value_threshold = tk.DoubleVar(value=0.05)
    log_fc_threshold = tk.DoubleVar(value=1.0)

    def select_file(var):
        file_path = filedialog.askopenfilename(
            title="Select File",
            filetypes=[("Excel Files", "*.xlsx"), ("All Files", "*.*")]
        )
        if file_path:
            var.set(file_path)

    def run_analysis():
        try:
            params = {
                "aging_file": aging_file.get(),
                "als_file": als_file.get(),
                "p_value_threshold": p_value_threshold.get(),
                "log_fc_threshold": log_fc_threshold.get(),
            }
            output_folder = filedialog.askdirectory(title="Select Output Folder")
            if not output_folder:
                raise ValueError("No output folder selected!")

            venn_path = analyze_data(params, output_folder)
            messagebox.showinfo("Success", f"Analysis completed! Venn diagram saved at {venn_path}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    # GUI Layout
    tk.Label(root, text="Load Aging Genes File:").pack(pady=5)
    tk.Entry(root, textvariable=aging_file, width=50).pack(pady=5)
    tk.Button(root, text="Browse", command=lambda: select_file(aging_file)).pack(pady=5)

    tk.Label(root, text="Load ALS Genes File:").pack(pady=5)
    tk.Entry(root, textvariable=als_file, width=50).pack(pady=5)
    tk.Button(root, text="Browse", command=lambda: select_file(als_file)).pack(pady=5)

    tk.Label(root, text="P-value Threshold:").pack(pady=5)
    tk.Entry(root, textvariable=p_value_threshold, width=10).pack(pady=5)

    tk.Label(root, text="Log2 Fold Change Threshold:").pack(pady=5)
    tk.Entry(root, textvariable=log_fc_threshold, width=10).pack(pady=5)

    tk.Button(root, text="Run Analysis", command=run_analysis).pack(pady=20)

    root.mainloop()

if __name__ == "__main__":
    start_gui()
