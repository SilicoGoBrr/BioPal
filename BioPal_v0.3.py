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

fasta_file = None
dir_path = None

Entrez.email = "" #INSTERT YOUR NCBI ACCOUNT MAIL (only neccessary for TaxaSage script)

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

def header_resumer():
    if not fasta_file:
        messagebox.showwarning("No File Selected", "Please select a FASTA file first using the 'Select Input File' button.")
        return

    sequences = SeqIO.parse(fasta_file, "fasta")
    new_headers = {}
    new_sequences = []
    header_map = []

    for record in sequences:
        full_header = record.description
        info = full_header.split("[organism=")[1].split("]")[0].split()
        new_header = info[0][0] + info[1][:5]  # Use at most 5 letters from the second word
        if new_header in new_headers:
            count = new_headers[new_header] + 1
            new_headers[new_header] = count
            new_header += str(count + 1)  # Append count + 1 to skip "1"
        else:
            new_headers[new_header] = 0
        new_record = SeqIO.SeqRecord(record.seq, id=new_header, description="")
        new_sequences.append(new_record)
        header_map.append((full_header, new_header))

    with open(os.path.join(dir_path, "header_data.csv"), "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Full header", "Resumed header"])
        for full_header, new_header in header_map:
            writer.writerow([full_header, new_header])

    with open(os.path.join(dir_path, "FASTA_con_headers_cortos.fasta"), "w") as handle:
        SeqIO.write(new_sequences, handle, "fasta")

def protparam_calculator():
    if not fasta_file:
        messagebox.showwarning("No File Selected", "Please select a FASTA file first using the 'Select Input File' button.")
        return

    with open(os.path.join(dir_path, "Propatam_Output.csv"), "w") as results:
        handle = open(fasta_file, "r")
        header = "Accession, Aa Num, MW (kDa), IP, GRAVY index, Pro, Asp, Lys, Asn, Glu, Arg, Sequence \n"
        results.write(header)

        for record in SeqIO.parse(handle, "fasta"):
            acc = str(record.id) + ","
            seq = str(record.seq.replace("X", ""))
            param = ProtParam.ProteinAnalysis(seq)
            num = str(len(seq)) + ","
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

def getFoldindex(seq):
    if not fasta_file:
        messagebox.showwarning("No File Selected", "Please select a FASTA file first using the 'Select Input File' button.")
        return
    
    r = requests.get("https://fold.proteopedia.org/cgi-bin/findex?m=json&sq=" + seq)
    data = json.loads(r.text)
    return data["findex"]

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

def get_taxonomic_info(organism):
    """
    Retrieves taxonomic information (Division, Order, Class, Family) for the given organism from NCBI Taxonomy.
    """
    try:
        search = Entrez.esearch(db="taxonomy", term=organism)
        result = Entrez.read(search)
        if not result["IdList"]:
            return None

        tax_id = result["IdList"][0]
        summary = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        tax_info = Entrez.read(summary)[0]

        division, order, class_, family = "N/A", "N/A", "N/A", "N/A"
        for lineage in tax_info.get("LineageEx", []):
            if lineage["Rank"] == "division":
                division = lineage["ScientificName"]
            elif lineage["Rank"] == "order":
                order = lineage["ScientificName"]
            elif lineage["Rank"] == "class":
                class_ = lineage["ScientificName"]
            elif lineage["Rank"] == "family":
                family = lineage["ScientificName"]

        return division, order, class_, family

    except Exception as e:
        return None

def phylogenetic_data_retrieval():
    """
    Retrieves Division, Order, Class, and Family information for each unique organism in the FASTA file.
    The output is written to a CSV file in the same directory as the input FASTA file.
    """
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
