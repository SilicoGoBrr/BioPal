import tkinter as tk
from tkinter import filedialog as fd, messagebox
import os
from Bio import SeqIO
from Bio import Entrez
from Bio.SeqUtils import ProtParam
import csv
import requests
import json
import time

# Path reset
fasta_file = None
dir_path = None

Entrez.email = "" #INSTERT YOUR NCBI ACCOUNT MAIL (only neccessary for TaxaSage script)

# Function for main input file selection (fasta, fa, faa, txt in fasta format)
def input_selector():
    global fasta_file, dir_path
    root.withdraw()  # Hide the root window while selecting file
    file_path = fd.askopenfilename(filetypes=[("FASTA files", "*.fasta *.fa")])
    if file_path:
        fasta_file = file_path
        dir_path = os.path.dirname(file_path)
        messagebox.showinfo("File Selected", f"Selected file: {fasta_file}")
    root.deiconify()  # Show the root window again after file selection

def split_fasta_file():
    if not fasta_file:
        messagebox.showwarning("No File Selected", "Please select a FASTA file first using the 'Select Input File' button.")
        return

    MAX_SEQS_PER_FILE = 100  # maximum number of sequences per output file
    with open(fasta_file, 'r') as f:
        seq_count = 0
        file_count = 1
        output_file = f'{os.path.splitext(fasta_file)[0]}_{file_count}.fasta'
        for line in f:
            if line.startswith('>'):
                # start of new sequence
                seq_count += 1
                if seq_count > MAX_SEQS_PER_FILE:
                    # start a new output file
                    file_count += 1
                    seq_count = 1
                    output_file = f'{os.path.splitext(fasta_file)[0]}_{file_count}.fasta'
            with open(output_file, 'a') as out_f:
                out_f.write(line)

# Function for header shortening.
def header_resumer():
    if not fasta_file:
        messagebox.showwarning("No File Selected", "Please select a FASTA file first using the 'Select Input File' button.")
        return
    # Set the FASTA file for parsing via SeqIO.
    sequences = SeqIO.parse(fasta_file, "fasta")
    new_headers = {}
    new_sequences = []
    header_map = []
    # For every sequence found in the file
    for record in sequences:
        #Store the header, find the "organism" tag, and extract the full extent of the species name.
        full_header = record.description
        info = full_header.split("[organism=")[1].split("]")[0].split()
        #Build a new resumed header using the first letter of the genus name and the first letters of the spp name.
        new_header = info[0][0] + info[1][:5]  # Use at most 5 letters from the second word
        #Check if such resumed header has been built before, and in that case add a number to the end to avoid recursive names. Otherwise leave it as is.
        if new_header in new_headers:
            count = new_headers[new_header] + 1
            new_headers[new_header] = count
            new_header += str(count + 1)  # Append count + 1 to skip "1"
        else:
            new_headers[new_header] = 0
        #Append the results to both the new resumed FASTA and the HEADER map csv file.
        new_record = SeqIO.SeqRecord(record.seq, id=new_header, description="")
        new_sequences.append(new_record)
        header_map.append((full_header, new_header))
    #Save the files in the FASTA file origin folder/path
    with open(os.path.join(dir_path, "header_data.csv"), "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Full header", "Resumed header"])
        for full_header, new_header in header_map:
            writer.writerow([full_header, new_header])

    with open(os.path.join(dir_path, "FASTA_con_headers_cortos.fasta"), "w") as handle:
        SeqIO.write(new_sequences, handle, "fasta")

# Function for protein parameter calculation
def protparam_calculator():
    if not fasta_file:
        messagebox.showwarning("No File Selected", "Please select a FASTA file first using the 'Select Input File' button.")
        return
    #Define output file
    with open(os.path.join(dir_path, "Propatam_Output.csv"), "w") as results:
        handle = open(fasta_file, "r")
        #Predefined (hard-coded) output header, If any modifications are wanted, you can find out how to add it in the Bio package (SeqUtils Protparam)
        header = "Accession, Aa Num, MW (kDa), IP, GRAVY index, Pro, Asp, Lys, Asn, Glu, Arg, Sequence \n"
        results.write(header)
        #Analyze and store results
        for record in SeqIO.parse(handle, "fasta"):
            acc = str(record.id) + ","
            seq = str(record.seq.replace("X", ""))
            param = ProtParam.ProteinAnalysis(seq)
            num = str(len(seq)) + ","
            #Express MW as kDa
            mw = str("%0.5f" % (param.molecular_weight()/1000)) +","
            ip = str(param.isoelectric_point()) +","
            gravy = str(param.gravy()) +","
            pro = str(param.count_amino_acids()['P']) + ","
            asp = str(param.count_amino_acids()['D']) + ","
            lys = str(param.count_amino_acids()['K']) + ","
            asn = str(param.count_amino_acids()['N']) + ","
            glu = str(param.count_amino_acids()['E']) + ","
            arg = str(param.count_amino_acids()['R']) + ","
            seqt = seq + "\n"
            all = acc + num + mw + ip + gravy + pro + asp + lys + asn + glu + arg + seqt
            results.write(all)

# Fold index calculator function. It uses INTERNET, as it makes requests.get calls from proteopedia's url. It's been deprecated by its non-online successor, but remains in code for reference and alternative use.
def getFoldindex(seq):
    if not fasta_file:
        messagebox.showwarning("No File Selected", "Please select a FASTA file first using the 'Select Input File' button.")
        return
    
    r = requests.get("https://fold.proteopedia.org/cgi-bin/findex?m=json&sq=" + seq)
    data = json.loads(r.text)
    return data["findex"]
# Successor to fold_index calculator function getFoldindex. Offline calculator.
def fold_index_calculator():
    if not fasta_file:
        messagebox.showwarning("No File Selected", "Please select a FASTA file first using the 'Select Input File' button.")
        return

    with open(fasta_file, "r") as fafa:
        for record in SeqIO.parse(fafa, "fasta"):
            sequence = str(record.seq)
            foldindex = getFoldindex(sequence)
            with open(os.path.join(dir_path, "FoldIndex_Output.csv"), "a") as results:
                results.write(str(record.id) + "," + str(foldindex) + '\n')

# This function requires internet as it fetches data from ENTREZ servers.
def phylogenetic_data_retrieval():
#  Retrieves Division, Order, Class, and Family information for each unique organism in the FASTA file.
#    The output is written to a CSV file in the same directory as the input FASTA file.
    if not fasta_file:
        messagebox.showwarning("No File Selected", "Please select a FASTA file first using the 'Select Input File' button.")
        return

    organisms = set()  # To avoid duplicates
    results = []

    try:
        with open(fasta_file, "r") as infile:
            sequences = SeqIO.parse(infile, "fasta")
            for record in sequences:
                header = record.description
                if "[organism=" in header:
                    organism = header.split("[organism=")[1].split("]")[0]
                    if organism not in organisms:
                        organisms.add(organism)
                        tax_info = get_taxonomic_info(organism)
                        if tax_info:
                            results.append([organism] + list(tax_info))
                        else:
                            results.append([organism, "N/A", "N/A", "N/A", "N/A"])
                        time.sleep(0.34)  # Throttle API requests to avoid rate-limiting

        # Write results to CSV
        output_file = os.path.join(dir_path, "phylogenetic_data.csv")
        with open(output_file, "w", newline="") as outfile:
            writer = csv.writer(outfile)
            writer.writerow(["Organism", "Division", "Order", "Class", "Family"])
            writer.writerows(results)

        messagebox.showinfo("Success", f"Phylogenetic data saved to {output_file}")

    except Exception as e:
        messagebox.showerror("Error", f"An error occurred while retrieving phylogenetic data: {str(e)}")

#Function fo fetch gff3 info from NCBI servers. Requires the input to include position in the header (this is available in Datasets NCBI)
def fetch_gff3(accession, start, stop, retries=3, timeout=30):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nuccore",
        "id": accession,
        "seq_start": start,
        "seq_stop": stop,
        "rettype": "gff",
        "retmode": "text"
    }

    for attempt in range(retries):
        try:
            response = requests.get(base_url, params=params, timeout=timeout)
            if response.status_code == 200:
                return response.text
            else:
                print(f"Error: Received status code {response.status_code}")
        except requests.exceptions.Timeout:
            print(f"Request timed out. Retrying ({attempt + 1}/{retries})...")
            time.sleep(5)
        except requests.exceptions.RequestException as e:
            print(f"An error occurred: {e}")
            break
    failed_accessions.append(accession)
    return None
# Parsing the gff3 information
def parse_gff3(gff3_data, organism, accession):
    genes = []
    current_gene = {
        'Species': organism,
        'GeneID': None,
        'Accession': None,
        'Description': None
    }
    capture_mrna = False

    for line in gff3_data.splitlines():
        if line.startswith("     gene"):
            # New gene section, reset current_gene
            current_gene = {
                'Species': organism,
                'Accession': None, # Reset Accession to None
                'GeneID': None,  # Reset GeneID to None
                'Description': None  # Reset Description to None
            }
            capture_mrna = False  # Reset capture flag for mRNA

        elif line.startswith("     mRNA"):
            # mRNA section - prepare to capture
            capture_mrna = True

        elif capture_mrna:
            # Look for the geneID and transcript_id lines within mRNA section
            if "/transcript_id=" in line:
                try:
                    transcript_id = line.split("/transcript_id=")[1].strip().strip('"')
                    print(transcript_id)
                    current_gene['Accession'] = transcript_id
                except IndexError:
                    current_gene['Accession'] = None  # If transcript_id is not found, set to None
            elif "/db_xref=" in line:
                try:
                    geneID = line.split('/db_xref="GeneID:')[1].strip().strip('"')
                    print(geneID)
                    current_gene['GeneID'] = geneID
                except IndexError:
                    current_gene['GeneID'] = None  # If GeneID is not found, set to None
            elif "/product=" in line:
                try:
                    description = line.split('/product="')[1].strip().strip('"')
                    current_gene['Description'] = description
                    print(description)
                except IndexError:
                    current_gene['Description'] = None  # If product description is not found, set to None

            # Once all fields (GeneID, Accession, Description) are captured, add the gene to the list
            if current_gene['GeneID']:
                genes.append(current_gene.copy())
                # Reset current_gene for the next gene
                current_gene = {
                    'Species': organism,
                    'Accession': None,
                    'GeneID': None,
                    'Description': None
                }
                capture_mrna = False  # Done capturing this mRNA section

    return genes

def parse_fasta_and_download_gff(fasta_file, output_csv):
    last_species = None
    with open(output_csv, mode='w', newline='') as csv_file:
        fieldnames = ['Species', 'Accession', 'Description', 'GeneID']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()

        for record in SeqIO.parse(fasta_file, "fasta"):
            header = record.description
            print(f"Header: {header}")

            try:
                # Extract accession, start, stop, and organism information from header
                parts = header.split()
                accession_region = parts[0]
                accession, region = accession_region.split(":")
                start, stop = region.lstrip("c").split("-")
                start, stop = int(start), int(stop)
                organism = header.split("[organism=")[1].split("]")[0]

                print(f"Parsed: Accession: {accession}, Start: {start}, Stop: {stop}, Organism: {organism}")
            except Exception as e:
                print(f"Error parsing header: {e}")
                continue

            if last_species and last_species != organism:
                # Write a blank line only if the species changes
                writer.writerow({})

            # Calculate upstream and downstream regions
            upstream_start = max(start - 20000, 0) # Change the 20000 integer to whatever number of bases upstream you want to analyze
            upstream_stop = start
            downstream_start = stop
            downstream_stop = stop + 20000 # Change the integer to whatever number of bases downstream you want to analyze

            # Fetch upstream GFF3 data
            gff3_data_upstream = fetch_gff3(accession, upstream_start, upstream_stop)
            if gff3_data_upstream:
                upstream_genes = parse_gff3(gff3_data_upstream, organism, accession)
                if upstream_genes:
                    for gene in upstream_genes:
                        writer.writerow(gene)
                else:
                    writer.writerow({'Species': organism, 'Accession': '', 'Description': 'No upstream genes found'})
            else:
                writer.writerow({'Species': organism, 'Accession': '', 'Description': 'No upstream GFF3 data retrieved'})

            # Write the query gene information after upstream genes
            query_gene = {
                'Species': organism,
                'Accession': accession,
                'Description': "Query Gene",
                'GeneID': '',  # Leave empty
            }
            writer.writerow(query_gene)

            # Fetch downstream GFF3 data
            gff3_data_downstream = fetch_gff3(accession, downstream_start, downstream_stop)
            if gff3_data_downstream:
                downstream_genes = parse_gff3(gff3_data_downstream, organism, accession)
                if downstream_genes:
                    for gene in downstream_genes:
                        writer.writerow(gene)
                else:
                    writer.writerow({'Species': organism, 'Accession': '', 'Description': 'No downstream genes found'})
            else:
                writer.writerow({'Species': organism, 'Accession': '', 'Description': 'No downstream GFF3 data retrieved'})

            last_species = organism
            writer.writerow({})

def select_fasta_file():
    root = tk.Tk()
    root.withdraw()
    file_path = fd.askopenfilename(
        title="Select the FASTA file",
        filetypes=[("FASTA files", "*.fasta"), ("FA files", "*.fa"), ("All files", "*.*")]
    )
    return file_path

def microsintenic_analyzer():
    fasta_file = select_fasta_file()

    if not fasta_file:
        print("No file selected.")
        return

    input_dir = os.path.dirname(fasta_file)
    output_csv = os.path.join(input_dir, "gff3_output.csv")

    print(f"Processing file: {fasta_file}")
    print(f"Output will be saved to: {output_csv}")

    parse_fasta_and_download_gff(fasta_file, output_csv)

    print("Processing complete.")

def show_help():
    messagebox.showinfo("Help", "This kit includes:\n"
                                 "1. Split FASTA file: This function takes any FASTA file provided and turns it into several files with a limit of 99 sequences in each file. Each file is numbered to keep track.\n"
                                 "2. Header resumer: This function grabs any FASTA file containing the species between brackets, like NCBI Database usually allows to download files, (e.g., XP_0000123.1 [organism=Ananas comosus]). The function then takes the first letter of the first word (A) and five from the second word (coreu) to build a new short name (Acoreu). It keeps track of repeated species by adding numbers after the first one (e.g., Acomos, Acomos2, Acomos3, etc.). It also returns a CSV file containing the original headers along with the new ones, so the user can decide to replace them again in the future if they want to.\n"
                                 "3. ProtParam calculator: This function takes the FASTA file and returns several parameters of interest from each sequence, as if bulk querying protparam Expasy tools. The output is a CSV file, which can be easily imported into Excel.\n"
                                 "Note: This program IGNORES 'X' characters in all sequences to perform calculations without errors.\n"
                                 "4. FoldIndex calculator: Queries proteopedia's fold tool and retrieves the fold index of the protein sequence.\n"
                                 "5. TaxaSage: Using the organism name taken from your selected fasta file (must include '[organism=' tag), uses NCBIs API Entrez services to gather Order, Class and Family and returns a .csv file for further use. \n"
                                 "6. Microsintenic retriever: Selecting a NCBI Database gene fasta file, it parses the data, collects gff3 data from NCBI servers and presents all coding genes surrounding the gene of interest both upstream and downstream in a readable CSV for the user. The distance is set to 20 kbp upstream from the gene start position, adn downstream from the end position."
                                 "\n"
                                 "BioPal v0.3.\n"
                                 "Disclaimer: This program is provided 'as is' without any guarantees or warranty. In connection with the program, I make no warranties of any kind, either express or implied, including but not limited to warranties of merchantability, fitness for a particular purpose, of title, or of noninfringement of third party rights. Use of the program is at your own risk, and I am not responsible for any misuse or the results produced by it. The program is subject to change at any moment without warning, and different versions of the program might produce slightly different results. Some scripts require an active internet connection to function correctly.")

def exit_func():
    if messagebox.askyesno("Exit", "Do you want to exit the program?"):
        root.destroy()

def main():
    global root
    root = tk.Tk()
    root.title("BioPal v0.3")

    frame = tk.Frame(root)
    frame.pack(padx=10, pady=10)

    buttons = [
        ("Select Input File", input_selector),
        ("Split FASTA File (Nuc/Prot)", split_fasta_file),
        ("Header Resumer (Nuc/Prot)", header_resumer),
        ("ProtParam Calculator (Prot)", protparam_calculator),
        ("Fold Index Calculator (Prot)", fold_index_calculator),
        ("Taxa Sage (Nuc/prot)", phylogenetic_data_retrieval),
        ("Microsintenic retriever (Nuc)", microsintenic_analyzer),
        ("Help", show_help),
        ("Exit", exit_func)
    ]

    for text, command in buttons:
        button = tk.Button(frame, text=text, command=command, width=30, height=2)
        button.pack(pady=5)

    root.mainloop()

if __name__ == "__main__":
    main()
