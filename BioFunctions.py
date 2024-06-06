import tkinter as tk
from tkinter import filedialog as fd
import os
from Bio import SeqIO
from Bio.SeqUtils import ProtParam
import csv
import requests
import json


def split_fasta_file():
    window = tk.Tk()
    file_path = fd.askopenfilename()

    MAX_SEQS_PER_FILE = 100  # maximum number of sequences per output file
    with open(file_path, 'r') as f:
        seq_count = 0
        file_count = 1
        output_file = f'{os.path.splitext(file_path)[0]}_{file_count}.fasta'
        for line in f:
            if line.startswith('>'):
                # start of new sequence
                seq_count += 1
                if seq_count > MAX_SEQS_PER_FILE:
                    # start a new output file
                    file_count += 1
                    seq_count = 1
                    output_file = f'{os.path.splitext(file_path)[0]}_{file_count}.fasta'
            with open(output_file, 'a') as out_f:
                out_f.write(line)
    return

def header_resumer():
    window = tk.Tk()
    file_path = fd.askopenfilename()
    # Get the directory of the input file
    dir_path = os.path.dirname(file_path)

    # Conversión del archivo FASTA a una lista
    sequences = SeqIO.parse(file_path, "fasta")

    # Create a dictionary to keep track of the new headers and another list to store the sequences
    new_headers = {}
    new_sequences = []

    # Iterate through the sequences in the FASTA file
    for record in sequences:
        # Get the full header
        full_header = record.description
        # Extract the relevant information: split the header, and use its information to create a 6 character long new header.
        info = full_header.split("[")[1].split("]")[0].split()
        new_header = info[0][0] + info[1][:6]
        # Check if the new header already exists, to see if it is necessary to add numeration.
        if new_header in new_headers:
            # Increment the count for this header, because it already exists
            count = new_headers[new_header] + 1
            new_headers[new_header] = count
            new_header += str(count)
        # The header does not exist, so just keep the header.
        else:
            new_headers[new_header] = 0
        # Create new SeqRecord with the same sequence data and new header
        new_record = SeqIO.SeqRecord(record.seq, id=new_header, description="")
        new_sequences.append(new_record)

    # Open a .csv file for writing
    with open(os.path.join(dir_path, "header_data.csv"), "w") as f:
        # Create a CSV writer
        writer = csv.writer(f)
        # Write the column headers
        writer.writerow(["Full header", "Resumed header"])
        # Write the header data to the .csv file
        for record in new_sequences:
            writer.writerow([full_header, record.id])

    # Write the modified sequences to a new FASTA file
    with open(os.path.join(dir_path, "FASTA con headers cortos.fasta"), "w") as handle:
        SeqIO.write(new_sequences, handle, "fasta")

    return

def protparam_calculator():
    #Open dialog for file selection
    window = tk.Tk()
    file_path = fd.askopenfilename()
    dir_path = os.path.dirname(file_path)
    #The optput file will contain your results in a CSV format with a proper header:
    with open(os.path.join(dir_path, "Propatam Output.csv"), "w") as results:
        #This will just call the FASTA file we need to make the bulk analysis:
        handle = open(file_path, "r")
        #This will be the headerfor our csv file:
        header = "Accession, Aa Num, MW (kDa), IP, GRAVY index, Pro, Asp, Lys, Asn, Glu, Arg, Sequence \n"
        results.write(header)

        for record in SeqIO.parse(handle, "fasta"):
            acc = str(record.id) + ","
            seq = str(record.seq)
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
        r = requests.get("https://fold.proteopedia.org/cgi-bin/findex?m=json&sq=" + seq)
        data = json.loads(r.text)
        return data["findex"]

def fold_index_calculator():
    window = tk.Tk()
    fasta_file = fd.askopenfilename()
    dir_path = os.path.dirname(fasta_file)
    with open(fasta_file, "r") as fafa:
        for record in SeqIO.parse(fafa, "fasta"):
            sequence = str(record.seq)
            foldindex = getFoldindex(sequence)
            with open(os.path.join(dir_path, "FoldIndex Output.csv"), "a") as results:
                results.write(str(record.id) + "," + str(foldindex) + '\n')

def input_selector():
    root = tk.Tk()
    fasta_file = fd.askopenfilename()
    dir_path = os.path.dirname(fasta_file)
    return fasta_file, dir_path

def main():
    root = tk.Tk()
    root.title("BioPal v0.1")

    frame = tk.Frame(root)
    frame.pack(padx=10, pady=10)
    select_input=tk.Button(frame, text="Select input file", command=input_selector)
    func1 = tk.Button(frame, text='Split FASTA File', command=split_fasta_file)
    func1.grid(row=3, column=0, padx=10, sticky="w")

    func2 = tk.Button(frame, text='Header resumer', command=header_resumer)
    func2.grid(row=4, column=0, padx=10, sticky="w")
    # Dictionary to map commands to functions
    commands = {
        '1': split_fasta_file(),
        '2': header_resumer(),
        '3': protparam_calculator(),
        '4': fold_index_calculator(),
        'help': help()
    }
    print('This tools are meant to be used with PROTEIN sequences, even though some of the functions might work with DNA sequences, reproducibility is not guaranteed.')
    while True:
        # Capture user input
        user_input = input("Enter command ('help' for commands): ").strip()
        
        if user_input.lower() == 'exit':
            print("Exiting program.")
            break
        elif user_input in commands:
            # Execute the corresponding function
            commands[user_input]()
        else:
            print("Unknown command. Write 'help' for available commands")
def help():
    print('This kit includes:')
    print('1. split_fasta_file(): This function takes any FASTA file provided and turns it in several files with a limit of 99 sequences in each file \n Each file is numbered to keep track.')
    print('2. header_resumer(): This functions grabs any fasta file containing the species between brackets, like NCBI ususally allows to download files, (e.g. XP_0000123.1 [Ananas comosus]) \n The function then takes the first letter of the first word (A) and five from the second word (coreu) to build a new short name (Acoreu)\n It keeps track of repeated species by adding numbers after the first one (e.g. Acoreu, Acoreu2, Acoreu3, etc.)\n It also returns a csv file containing the original headers along with the new ones, so the user can decide to replace them again in the future if he/she wants to.')
    print('3. protparam_calculator(): This function takes the FASTA file and return several parameter of interest from each sequence, as if bulk querying protparam Expasy tools. The output is a csv file, which can be easilty imported in excel.')
    print('NOTE: protparam_calculator() will not work if your sequences contain non-coding characters in their sequences (e.g. "X"). Be sure that you have removed orreplaced such characters.')
    print('4. fold_index_calculator(): Queryies proteopedias fold tool and retrieves the fold index of the protein sequence.')
    main()

if __name__ == "__main__":
    main()