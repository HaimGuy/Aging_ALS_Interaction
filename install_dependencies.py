import os
import subprocess
import sys

# List of required libraries
libraries = [
    "pandas",
    "numpy",
    "matplotlib",
    "seaborn",
    "scipy",
    "openpyxl",
    "tkintertable",  # For GUI, if required
    "rpy2",          # For R-Python integration
    "matplotlib-venn",  # For Venn diagrams
]

# Function to install a library
def install_library(library):
    try:
        print(f"Installing {library}...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", library])
        print(f"{library} installed successfully.")
    except subprocess.CalledProcessError:
        print(f"Error installing {library}. Please install it manually.")

# Function to check and install all libraries
def install_libraries():
    for library in libraries:
        install_library(library)

# Function to install external tools (optional)
def install_external_tools():
    # Example for installing system-wide tools like 'featureCounts'
    featurecounts_command = "conda install -c bioconda subread"
    try:
        print(f"Installing featureCounts using: {featurecounts_command}")
        subprocess.run(featurecounts_command, shell=True, check=True)
        print("featureCounts installed successfully.")
    except subprocess.CalledProcessError:
        print("Error installing featureCounts. Please install it manually.")

if __name__ == "__main__":
    print("Starting dependency installation...")
    install_libraries()
    print("Checking external dependencies...")
    install_external_tools()
    print("All dependencies installed successfully.")
