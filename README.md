# BioPal

**BioPal (v0.2)** is a bioinformatics toolkit designed to process FASTA sequence files. The tool provides several functionalities such as splitting FASTA files, calculating protein parameters, querying taxonomic information from NCBI, and much more. It uses the `tkinter` library to provide a user-friendly graphical interface for easy file input and function selection.

## Features

1. **Split FASTA File**: Splits a FASTA file into multiple smaller files with a maximum of 99 sequences per file. Sometimes required for post prrocessing.
   
2. **Header Resumer**: Resumes long headers into shorter, standardized ones (e.g., based on the organism name from the NCBI format `[organism=...]`) and outputs a CSV mapping the original and new headers. It provides both a new FASTA file with the new short headers and the sequences, and a CSV file with both the "old" and "new/short" header names for adequate tracking of sequences.

3. **ProtParam Calculator**: Performs bulk calculations of various protein properties (e.g., molecular weight, isoelectric point, etc.) similar to ExPASy's ProtParam tool and outputs the results into a CSV file. Note: This program IGNORES "X" characters in all sequences to perform calculations without errors. Returns a CSV file with the results. So far, this feature is still hard-coded, and the user can't change the output of the program.

4. **Fold Index Calculator**: Queries the proteopedia fold index tool for each sequence in the FASTA file and outputs the fold index of each sequence to a CSV file.

5. **Taxa Sage**: Queries taxonomic information (Division, Order, Class, Family) for organisms in the FASTA file (requires the presence of `[organism=...]` in the header) and writes the results to a CSV file.

6. **Help Menu**: Provides a description of the tool's functionalities.

7. **Exit**: Safely closes the application. Program does not hold path/file information.

## Requirements

This tool requires the following Python libraries to be installed:

- `tkinter` for the graphical user interface.
- `biopython` for parsing FASTA files and retrieving taxonomic information.
- `requests` and `json` for querying online databases like proteopedia and NCBI.

You can install the necessary dependencies using:

```bash
pip install biopython requests
```

## Installation

1. Clone the repository:

```bash
git clone https://github.com/SilicoGoBrr/BioPal.git
cd BioPal
```

2. Install the dependencies:

```bash
pip install -r requirements.txt
```

3. Run the program:

```bash
python biopal.py
```

## How to Use

1. **Select Input File**: Click the "Select Input File" button to choose your FASTA file.

2. **Choose an Operation**:
   - Split the FASTA file.
   - Resume headers.
   - Calculate protein parameters.
   - Calculate fold index.
   - Retrieve taxonomic information (Taxa Sage).

3. The results will be saved in the same directory as the input file, with appropriate file names based on the operation performed.

## Taxa Sage: Entrez API Configuration

The `Taxa Sage` feature uses NCBI’s Entrez API to retrieve taxonomic data. For this, you need to specify your email address, as required by NCBI's Entrez API.

In the code, locate the following line:

```python
Entrez.email = ""  # Add your email here
```

Replace it with your valid email address:

```python
Entrez.email = "johndoe@fakename.com"
```

This step is necessary for the Entrez API requests to work properly.

## Known Issues

- **Rate-Limiting**: The NCBI Entrez API may impose rate limits. To avoid being rate-limited, the tool introduces a short delay between API requests when using the Taxa Sage feature.

- **FASTA Format Requirements**: The input FASTA file must contain `[organism=...]` tags for the Taxa Sage function to work correctly.

## Disclaimer

This program is provided "as is" without any warranties or guarantees. Use it at your own risk. Some features may require an active internet connection.
