import tkinter as tk
from tkinter import filedialog as fd, messagebox
import pandas as pd
import os
from Bio import SeqIO
from Bio import Entrez
from Bio.SeqUtils import ProtParam
import csv
import requests
import json
import time
import matplotlib.pyplot as plt

# Path reset
fasta_file = None
dir_path = None

Entrez.email = "" #INSTERT YOUR NCBI ACCOUNT MAIL (only neccessary for TaxaSage script)

# Function for main input file selection (fasta, fa, faa, txt in fasta format)
def input_selector():
    global fasta_file, dir_path
    root.withdraw()  # Hide the root window while selecting file
    file_path = fd.askopenfilename(filetypes=[("FASTA files", "*.fasta *.fa"), ("Text files", "*.txt"), ("All files", "*.*")])
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
        header = "Accession, Aa Num, MW (kDa), IP, GRAVY index, Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His, Ile, Leu, Lys, Met, Phe, Pro, Ser, Thr, Trp, Tyr, Val, Sequence \n"
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
            ala = str(param.count_amino_acids()['A']) + ","
            arg = str(param.count_amino_acids()['R']) + ","
            asn = str(param.count_amino_acids()['N']) + ","
            asp = str(param.count_amino_acids()['D']) + ","
            cys = str(param.count_amino_acids()['C']) + ","
            gln = str(param.count_amino_acids()['Q']) + ","
            glu = str(param.count_amino_acids()['E']) + ","
            gly = str(param.count_amino_acids()['G']) + ","
            his = str(param.count_amino_acids()['H']) + ","
            ile = str(param.count_amino_acids()['I']) + ","
            leu = str(param.count_amino_acids()['L']) + ","
            lys = str(param.count_amino_acids()['K']) + ","
            met = str(param.count_amino_acids()['M']) + ","
            phe = str(param.count_amino_acids()['F']) + ","
            pro = str(param.count_amino_acids()['P']) + ","
            ser = str(param.count_amino_acids()['S']) + ","
            thr = str(param.count_amino_acids()['T']) + ","
            trp = str(param.count_amino_acids()['W']) + ","
            tyr = str(param.count_amino_acids()['Y']) + ","
            val = str(param.count_amino_acids()['V']) + ","
            seqt = seq + "\n"
            all = acc + num + mw + ip + gravy + ala + arg + asn + asp + cys + gln + glu + gly + his + ile + leu + lys + met + phe + pro + ser + thr + trp + tyr + val + seqt
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
            writer.writerow(["Organism", "Class", "Order", "Family"])
            writer.writerows(results)

        messagebox.showinfo("Success", f"Phylogenetic data saved to {output_file}")

    except Exception as e:
        messagebox.showerror("Error", f"An error occurred while retrieving phylogenetic data: {str(e)}")

def get_taxonomic_info(organism):
    """
    Retrieves taxonomic information (Class, Order, Family) for the given organism from NCBI Taxonomy.
    """
    try:
        search = Entrez.esearch(db="taxonomy", term=organism)
        result = Entrez.read(search)
        if not result["IdList"]:
            return None

        tax_id = result["IdList"][0]
        summary = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        tax_info = Entrez.read(summary)[0]

        order, class_, family = "N/A", "N/A", "N/A"
        for lineage in tax_info.get("LineageEx", []):
            if lineage["Rank"] == "class":
                order = lineage["ScientificName"]
            elif lineage["Rank"] == "order":
                class_ = lineage["ScientificName"]
            elif lineage["Rank"] == "family":
                family = lineage["ScientificName"]

        return order, class_, family

    except Exception as e:
        return None

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
    failed_accessions = []
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

def parse_fasta_and_download_gff(fasta_file, output_csv):
    last_species = None
    seen_genes = set()  # Set to track unique genes
    query_gene_info = None

    with open(output_csv, mode='w', newline='') as csv_file:
        fieldnames = ['Species', 'Accession', 'Description', 'GeneID', 'Start', 'Stop']
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
            upstream_start = max(start - 30000, 0)  # Default 30000 upstream
            upstream_stop = start
            downstream_start = stop + 1
            downstream_stop = stop + 30000  # Default 30000 downstream

            # Fetch upstream GFF3 data
            gff3_data_upstream = fetch_gff3(accession, upstream_start, upstream_stop)
            if gff3_data_upstream:
                upstream_genes = parse_gff3(gff3_data_upstream, organism, accession)
                if upstream_genes:
                    for gene in upstream_genes:
                        gene_tuple = (gene['Accession'], gene['Description'], gene['GeneID'])
                        if gene_tuple not in seen_genes:
                            writer.writerow(gene)
                            seen_genes.add(gene_tuple)
                else:
                    writer.writerow({'Species': organism, 'Accession': '', 'Description': 'No upstream genes found'})
            else:
                writer.writerow({'Species': organism, 'Accession': '', 'Description': 'No upstream GFF3 data retrieved'})

            # Fetch downstream GFF3 data
            gff3_data_downstream = fetch_gff3(accession, downstream_start, downstream_stop)
            if gff3_data_downstream:
                downstream_genes = parse_gff3(gff3_data_downstream, organism, accession)

                if downstream_genes:
                    # Store the query gene information and write it once
                    query_gene_info = downstream_genes[0]
                    query_gene_tuple = (query_gene_info['Accession'], query_gene_info['Description'], query_gene_info['GeneID'])

                    if query_gene_tuple not in seen_genes:
                        writer.writerow(query_gene_info)
                        seen_genes.add(query_gene_tuple)

                    # Skip the first downstream gene (query gene) and write the rest
                    for gene in downstream_genes[1:]:
                        gene_tuple = (gene['Accession'], gene['Description'], gene['GeneID'])
                        if gene_tuple not in seen_genes:
                            writer.writerow(gene)
                            seen_genes.add(gene_tuple)

                    # If the combined gene list is less than or equal to 6 genes, extend the search region
                    if len(upstream_genes) + len(downstream_genes) <= 6:
                        print(f"Few genes found, expanding search for {accession}...")
                        upstream_start = max(start - 40000, 0)  # Expand upstream by an additional 10000 bases
                        downstream_stop = stop + 40000  # Expand downstream by 10000 bases
                        
                        # Fetch extended upstream GFF3 data
                        gff3_data_upstream = fetch_gff3(accession, upstream_start, upstream_stop)
                        if gff3_data_upstream:
                            upstream_genes = parse_gff3(gff3_data_upstream, organism, accession)
                            if upstream_genes:
                                for gene in upstream_genes:
                                    gene_tuple = (gene['Accession'], gene['Description'], gene['GeneID'])
                                    if gene_tuple not in seen_genes:
                                        writer.writerow(gene)
                                        seen_genes.add(gene_tuple)
                        
                        # Fetch extended downstream GFF3 data
                        gff3_data_downstream = fetch_gff3(accession, downstream_start, downstream_stop)
                        if gff3_data_downstream:
                            downstream_genes = parse_gff3(gff3_data_downstream, organism, accession)
                            if downstream_genes:
                                for gene in downstream_genes[1:]:  # Skip the first gene again in the extended query
                                    gene_tuple = (gene['Accession'], gene['Description'], gene['GeneID'])
                                    if gene_tuple not in seen_genes:
                                        writer.writerow(gene)
                                        seen_genes.add(gene_tuple)
                else:
                    writer.writerow({'Species': organism, 'Accession': '', 'Description': 'No downstream genes found'})
            else:
                writer.writerow({'Species': organism, 'Accession': '', 'Description': 'No downstream GFF3 data retrieved'})

            last_species = organism

# Parsing the gff3 information 
def parse_gff3(gff3_data, organism, accession):
    genes = []
    current_gene = {
        'Species': organism,
        'GeneID': None,
        'Accession': None,
        'Description': None,
        'Start': None,
        'Stop': None
    }
    capture_mrna = False

    for line in gff3_data.splitlines():
        if line.startswith("     gene"):
            # New gene section, reset current_gene
            current_gene = {
                'Species': organism,
                'Accession': None,  # Reset Accession to None
                'GeneID': None,  # Reset GeneID to None
                'Description': None,  # Reset Description to None
                'Start': None,  # Reset Start
                'Stop': None  # Reset Stop
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
                    current_gene['Accession'] = transcript_id
                except IndexError:
                    current_gene['Accession'] = None  # If transcript_id is not found, set to None
            elif "/db_xref=" in line:
                try:
                    geneID = line.split('/db_xref="GeneID:')[1].strip().strip('"')
                    current_gene['GeneID'] = geneID
                except IndexError:
                    current_gene['GeneID'] = None  # If GeneID is not found, set to None
            elif "/product=" in line:
                try:
                    description = line.split('/product="')[1].strip().strip('"')
                    current_gene['Description'] = description
                except IndexError:
                    current_gene['Description'] = None  # If product description is not found, set to None
            elif "/location=" in line:
                try:
                    start_stop = line.split('/location="')[1].strip().strip('"').split("..")
                    current_gene['Start'], current_gene['Stop'] = int(start_stop[0]), int(start_stop[1])
                except (IndexError, ValueError):
                    current_gene['Start'], current_gene['Stop'] = None, None  # If location is not found or not valid, set to None

            # Once all fields (GeneID, Accession, Description, Start, Stop) are captured, add the gene to the list
            if current_gene['GeneID']:
                genes.append(current_gene.copy())
                # Reset current_gene for the next gene
                current_gene = {
                    'Species': organism,
                    'Accession': None,
                    'GeneID': None,
                    'Description': None,
                    'Start': None,
                    'Stop': None
                }
                capture_mrna = False  # Done capturing this mRNA section

    # Sort the genes by the Start position before returning
    genes = sorted(genes, key=lambda gene: gene['Start'] if gene['Start'] is not None else float('inf'))
    return genes

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

def identify_blocks(df):
    """
    Identifies blocks of data in the DataFrame based on consecutive non-empty rows.
    """
    blocks = []
    current_block = []
    
    for _, row in df.iterrows():
        if not row.isnull().all():  # Check if the row is not entirely empty
            current_block.append(row)
        else:
            if current_block:
                blocks.append(pd.DataFrame(current_block))
                current_block = []
    
    if current_block:
        blocks.append(pd.DataFrame(current_block))
    
    return blocks

def compare_blocks(blocks, min_similar_genes, warning_file):
    """
    Compares gene blocks and groups similar ones, generating a warning for lower than required matches.
    """
    groups = []
    assigned = set()
    warnings = []  # Lista para almacenar las advertencias
    
    for i, block1 in enumerate(blocks):
        if i in assigned:
            continue
        
        current_group = [block1.copy()]
        assigned.add(i)
        matches = {}  # Diccionario para almacenar especies coincidentes por fila
        
        for j, block2 in enumerate(blocks):
            if i == j or j in assigned:
                continue
            
            common_descriptions = set(block1['Description'].unique()) & set(block2['Description'].unique())
            num_common = len(common_descriptions)
            
            # Si la cantidad de coincidencias supera el umbral mínimo, se agrupan los bloques
            if num_common >= min_similar_genes:
                current_group.append(block2.copy())
                assigned.add(j)
                species_name = block2.iloc[0, 0]  # Obtener el nombre de la especie del bloque 2
                
                for index, row in block1.iterrows():
                    if row['Description'] in common_descriptions:
                        if index not in matches:
                            matches[index] = set()
                        matches[index].add(species_name)
            # Si las coincidencias son menores que el mínimo, se genera una advertencia
            elif num_common > 0:
                warnings.append([
                    f"Bloque {i+1}",
                    block1.iloc[0, 0],  # Especie del bloque 1
                    num_common,
                    ", ".join(common_descriptions),
                    block2.iloc[0, 0]  # Especie del bloque 2
                ])
        
        if len(current_group) > 1:  # Solo si el grupo tiene más de un bloque
            for index, species_set in matches.items():
                if 'Matching Species' not in block1.columns:
                    block1['Matching Species'] = None  # Inicializar columna si no existe
                block1.at[index, 'Matching Species'] = ', '.join(species_set)
            groups.append(current_group)
    
    # Escribir el archivo de advertencias si se encontraron coincidencias menores
    if warnings:
        with open(warning_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow([
                "Las siguientes especies mostraron coincidencias en menor número al indicado por el usuario:"
            ])
            writer.writerow(["Bloque", "Especie del bloque", "Número de coincidencias", "Descripción(s) que coinciden", "Especie del bloque con coincidencia"])
            for warning in warnings:
                writer.writerow(warning)
    
    return groups

def create_output_files(groups, file_format, output_dir):
    """
    Creates output files for each group in the selected format (CSV or XLSX), adding a counter for the number of blocks.
    """
    file_count = 0
    
    for idx, group in enumerate(groups):
        filename = os.path.join(output_dir, f"group_{idx + 1}.{file_format}")  # Asegura que el archivo se guarda en la misma carpeta
        
        if file_format == "csv":
            with open(filename, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(["Number of Blocks in Group:", len(group)])  # Add block count
                
                for block in group:
                    writer.writerow([])  # Add empty line for readability
                    writer.writerow(block.columns)
                    for _, row in block.iterrows():
                        writer.writerow(row.values)
        
        elif file_format == "xlsx":
            with pd.ExcelWriter(filename) as writer:
                for b_idx, block in enumerate(group):
                    sheet_name = f"Block_{b_idx+1}"
                    block.to_excel(writer, sheet_name=sheet_name, index=False)
                writer.sheets['Sheet1'].cell(row=1, column=1, value=f"Number of Blocks in Group: {len(group)}")
        
        file_count += 1
    
    return file_count

def browse_file():
    global gfffile
    gfffile = fd.askopenfilename(filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
    entry.delete(0, 'end')
    entry.insert(0, gfffile)
    return gfffile

def start_analysis():
    global gfffile
    min_similar_genes = int(entry_min_genes.get())
    file_format = file_format_var.get()
    
    try:
        df = pd.read_csv(gfffile)
        blocks = identify_blocks(df)
        
        # Nombre del archivo de advertencias
        warning_file = os.path.join(os.path.dirname(gfffile), "Other matches.csv")
        
        groups = compare_blocks(blocks, min_similar_genes, warning_file)
        
        # Obtener el directorio de salida (misma carpeta que el archivo de entrada)
        output_dir = os.path.dirname(gfffile)  # Guardar archivos en el mismo directorio que el archivo de entrada
        print(f"Output directory: {output_dir}")  # Verifica la ruta

        file_count = create_output_files(groups, file_format, output_dir)
        
        messagebox.showinfo("Analysis Completed", f"Analysis finished. {file_count} files created. Warnings saved to {warning_file}")
    except FileNotFoundError:
        messagebox.showerror("Error", "File not found.")

# def create_similarity_graph(blocks, min_similar_genes):
#     """
#     Creates a graph where nodes represent species and edges represent similarity based on matching genes.
#     The thickness of edges depends on the number of matches, and the distance between nodes depends on similarity.
#     """
#     G = nx.Graph()  # Create an undirected graph
#     warnings = []  # List to store warnings if needed
    
#     for i, block1 in enumerate(blocks):
#         species1 = block1.iloc[0, 0]  # Get species name from the first column
#         G.add_node(species1)  # Add the node for this species
        
#         for j, block2 in enumerate(blocks):
#             if i >= j:  # To avoid duplicate comparisons
#                 continue
            
#             species2 = block2.iloc[0, 0]
#             common_descriptions = set(block1['Description'].unique()) & set(block2['Description'].unique())
#             num_common = len(common_descriptions)
            
#             if num_common >= min_similar_genes:
#                 # Add an edge between the nodes with the weight being the number of common descriptions
#                 G.add_edge(species1, species2, weight=num_common)
#             elif num_common > 0:
#                 warnings.append([
#                     f"Bloque {i+1}",
#                     species1,
#                     num_common,
#                     ", ".join(common_descriptions),
#                     species2
#                 ])
    
#     return G, warnings

# def plot_similarity_graph(G):
#     """
#     Plots the similarity graph using NetworkX and Matplotlib.
#     """
#     # Usar un layout que no agrupe demasiado los nodos
#     pos = nx.spring_layout(G, k=0.5, seed=42)  # Ajusta 'k' para más separación
    
#     # Ajuste de grosor de las líneas
#     edge_widths = [d['weight'] * 0.5 for u, v, d in G.edges(data=True)]
    
#     # Crear la figura
#     plt.figure(figsize=(14, 14))  # Aumentar el tamaño para mejor visibilidad
    
#     # Dibuja el gráfico con una paleta de colores ajustada
#     nx.draw(G, pos, with_labels=True, node_size=1000, node_color='skyblue', font_size=10,
#             font_weight='bold', width=edge_widths, edge_color=edge_widths, cmap=plt.cm.Blues)
    
#     plt.title("Gene Block Similarity Map", fontsize=16)
#     plt.show()

def ortho_finder_app():
    root = tk.Tk()
    root.title("Gene Block Comparator")

    label = tk.Label(root, text="Select the CSV file:")
    label.pack()
    global entry
    entry = tk.Entry(root, width=50)
    entry.pack()
    button = tk.Button(root, text="Select CSV File", command=browse_file)
    button.pack()

    label_min_genes = tk.Label(root, text="Minimum number of gene matches:")
    label_min_genes.pack()
    global entry_min_genes
    entry_min_genes = tk.Entry(root, width=10)
    entry_min_genes.pack()

    label_format = tk.Label(root, text="Select output format:")
    label_format.pack()
    global file_format_var
    file_format_var = tk.StringVar(value="csv")
    csv_button = tk.Radiobutton(root, text="CSV", variable=file_format_var, value="csv")
    csv_button.pack()
    xlsx_button = tk.Radiobutton(root, text="XLSX", variable=file_format_var, value="xlsx")
    xlsx_button.pack()

    button_start = tk.Button(root, text="Start Analysis", command=start_analysis)
    button_start.pack()

def show_help():
    messagebox.showinfo("Help", "This kit includes:\n"
                                 "1. Split FASTA file: This function takes any FASTA file provided and turns it into several files with a limit of 99 sequences in each file. Each file is numbered to keep track.\n"
                                 "2. Header resumer: This function grabs any FASTA file containing the species between brackets, like NCBI Database usually allows to download files, (e.g., XP_0000123.1 [organism=Ananas comosus]). The function then takes the first letter of the first word (A) and five from the second word (comos) to build a new short name (Acomos). It keeps track of repeated species by adding numbers after the first one (e.g., Acomos, Acomos2, Acomos3, etc.). It also returns a CSV file containing the original headers along with the new ones, so the user can decide to replace them again in the future if they want to.\n"
                                 "3. ProtParam calculator: This function takes the protein FASTA file and returns several parameters of interest from each sequence, as if bulk querying protparam Expasy tools. The output is a CSV file, which can be easily imported into Excel.\n"
                                 "Note: This program IGNORES 'X' characters in all sequences to perform calculations without errors.\n"
                                 "4. FoldIndex calculator: Queries proteopedia's FOldIndex tool and retrieves the fold index of the protein sequence. Note: This step is slow, as it queries proteins one-by-one\n"
                                 "5. TaxaSage: Using the organism name taken from your selected fasta file (must include '[organism=' tag), uses NCBIs API Entrez services to gather Order, Class and Family and returns a .csv file for further use. \n"
                                 "6. Microsintenic retriever: Selecting a NCBI Database gene fasta file, it parses the data, collects gff3 data from NCBI servers and presents all coding genes surrounding the gene of interest both upstream and downstream in a readable CSV for the user. The distance is set to 30 kbp upstream from the gene start position, and downstream from the end position. If few genes are detected, program automatically extends this by 10 kbps."
                                 "\n"
                                 "BioPal v0.5.\n"
                                 "Disclaimer: This program is provided 'as is' without any guarantees or warranty. In connection with the program, I make no warranties of any kind, either express or implied, including but not limited to warranties of merchantability, fitness for a particular purpose, of title, or of noninfringement of third party rights. Use of the program is at your own risk, and I am not responsible for any misuse or the results produced by it. The program is subject to change at any moment without warning, and different versions of the program might produce slightly different results. Some scripts require an active internet connection to function correctly.")

def exit_func():
    if messagebox.askyesno("Exit", "Do you want to exit the program?"):
        root.destroy()

def main():
    global root
    root = tk.Tk()
    root.title("BioPal v0.5")

    frame = tk.Frame(root)
    frame.pack(padx=5, pady=5)

    buttons = [
        ("Select Input File", input_selector),
        ("Split FASTA File", split_fasta_file),
        ("Header Resumer", header_resumer),
        ("ProtParam Calculator", protparam_calculator),
        ("Fold Index Calculator", fold_index_calculator),
        ("Taxa Sage", phylogenetic_data_retrieval),
        ("Microsyntenic retriever (MR) (Nuc)", microsintenic_analyzer),
        ("Ortho Finder (MR output)", ortho_finder_app),
        ("Help", show_help),
        ("Exit", exit_func)
    ]

    for text, command in buttons:
        button = tk.Button(frame, text=text, command=command, width=30, height=2)
        button.pack(pady=5)

    root.mainloop()

if __name__ == "__main__":
    main()
