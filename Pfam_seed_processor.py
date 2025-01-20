"""
Este script toma los alineamientos que se pueden obtener de las accesiones de PFam (dentro de InterPro) 
con el objetivo de "limpiarlas" y convertirlas en un FASTA conteniendo las secuencias sin alinear.
El objetivo es poder utilizarlas para realizar BLASTs y personalizar los parámetros de búsqueda en el mismo
para aprofundizar el análisis en secuencias similares.
"""

import tkinter as tk
from tkinter import filedialog as fd, messagebox
import os

# Path reset
seed_file = None
dir_path = None

# Selector de archivo seed
def input_selector():
    global seed_file, dir_path
    file_path = fd.askopenfilename(filetypes=[("Seed FILES", "*.seed"), ("Text files", "*.txt"), ("All files", "*.*")])
    if file_path:
        seed_file = file_path
        dir_path = os.path.dirname(seed_file)
        messagebox.showinfo("File Selected", f"Selected file: {seed_file}")
    return seed_file

# Procesamiento del archivo
def seed_processor(file_path):
    fasta_resultado = []
    try:
        with open(file_path, "r") as file:
            for linea in file:
                # Ignorar líneas que comiencen con '#'
                if linea.startswith("#"):
                    continue
                # Separar el identificador y la secuencia
                partes = linea.split()
                if len(partes) == 2:
                    identificador, secuencia = partes
                    # Limpiar caracteres no alfabéticos de la secuencia
                    secuencia_limpia = ''.join(filter(str.isalpha, secuencia))
                    # Agregar el formato FASTA a la lista
                    fasta_resultado.append(f">{identificador}")
                    fasta_resultado.append(secuencia_limpia)
    except Exception as e:
        messagebox.showerror("Error", f"Error al procesar el archivo: {str(e)}")
        return []

    # Imprimir todas las secuencias en formato FASTA al final
    if fasta_resultado:
        print("\n".join(fasta_resultado))
    else:
        print("No se encontraron secuencias válidas.")
    
    return fasta_resultado

# Main
if __name__ == "__main__":
    input_selector()  # Seleccionar el archivo
    if seed_file:
        seed_processor(seed_file)  # Procesar el archivo seleccionado
