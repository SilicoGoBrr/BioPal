import csv
import requests
from Bio import Entrez

# Configure Entrez email
Entrez.email = "axeljrizzo@gmail.com"  # Replace with your email

def get_taxonomic_info(organism):
    """
    Retrieves taxonomic information (Division, Order, Class, Family) for the given organism from NCBI Taxonomy.
    """
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

def process_input_file(input_file, output_file):
    """
    Processes the input file to extract organism names and retrieve taxonomic data for each unique organism.
    Writes the results to an output CSV file.
    """
    organisms = set()  # To avoid duplicates
    results = []

    with open(input_file, "r") as infile:
        reader = csv.reader(infile)
        next(reader)  # Skip the header
        for row in reader:
            full_header = row[0]
            if "[organism=" in full_header:
                organism = full_header.split("[organism=")[1].split("]")[0]
                if organism not in organisms:
                    organisms.add(organism)
                    tax_info = get_taxonomic_info(organism)
                    if tax_info:
                        results.append([organism] + list(tax_info))
                    else:
                        results.append([organism, "N/A", "N/A", "N/A", "N/A"])

    # Write results to CSV
    with open(output_file, "w", newline="") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(["Organism", "Division", "Order", "Class", "Family"])
        writer.writerows(results)

    print(f"Taxonomic data written to {output_file}")

if __name__ == "__main__":
    input_file = "input_file.csv"  # Replace with your input file path
    output_file = "phylogenetic_data.csv"  # Replace with your desired output file path

    process_input_file(input_file, output_file)
