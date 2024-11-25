import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

# Diccionario del código genético
CODON_TABLE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}

# Función para traducir una secuencia de ADN a proteínas
def translate_dna_to_protein(dna_sequence):
    proteins = []
    sequence_length = len(dna_sequence)
    for i in range(0, sequence_length, 3):  # Recorre la secuencia en pasos de 3 (longitud de un codón)
        codon = dna_sequence[i:i+3]
        if len(codon) == 3:  # Verifica que sea un codón completo
            protein = CODON_TABLE.get(codon, '')
            if protein == '_':  # Codón de parada
                break
            proteins.append(protein)
    return ''.join(proteins)

# Función para contar las proteínas
def count_proteins(protein_sequence):
    return Counter(protein_sequence)

# Función para contar los nucleótidos
def count_nucleotides(dna_sequence):
    return Counter(dna_sequence)

# Función para crear un gráfico de barras proteínas
def plot_protein_counts(protein_counts, title):
    fig, ax = plt.subplots()
    proteins = list(protein_counts.keys())
    counts = list(protein_counts.values())
    
    ax.bar(proteins, counts, color='skyblue')
    ax.set_xlabel("Tipo de Proteína")
    ax.set_ylabel("Cantidad")
    ax.set_title(title)
    return fig

# Función para crear un gráfico de barras nucleótidos
def plot_counts(counts, title, xlabel, ylabel, color='skyred'):
    fig, ax = plt.subplots()
    labels = list(counts.keys())
    values = list(counts.values())
    
    ax.bar(labels, values, color=color)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    return fig

# Función para crear un gráfico de barras apiladas
def plot_stacked_bar(data1, data2, labels, title, xlabel, ylabel, color1='skyblue', color2='orange'):
    fig, ax = plt.subplots()
    indices = np.arange(len(labels))
    
    # Obtener valores
    values1 = [data1.get(label, 0) for label in labels]
    values2 = [data2.get(label, 0) for label in labels]
    
    # Gráficas apiladas
    ax.bar(indices, values1, color=color1, label='Secuencia 1')
    ax.bar(indices, values2, bottom=values1, color=color2, label='Secuencia 2')
    
    ax.set_xticks(indices)
    ax.set_xticklabels(labels)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    return fig

# Configuración de Streamlit
st.title("Dashboard de Traducción de ADN a Proteínas")
st.markdown("""
Este dashboard permite introducir una secuencia de nucleótidos (ADN) para traducirla 
en proteínas y contar la cantidad de cada una presente en la secuencia.
""")

st.sidebar.title("Conceptos")
st.sidebar.write("""
**Proteínas**: Las proteínas son moléculas formadas por cadenas de aminoácidos, que a su vez son codificadas por las secuencias de ADN. Estas proteínas realizan una variedad de funciones dentro de las células, como enzimas, estructuras y señales.

**Nucleótidos**: Los nucleótidos son las unidades básicas del ADN. Cada nucleótido está formado por un azúcar, un fosfato y una base nitrogenada (adenina, timina, citosina o guanina). Las secuencias de nucleótidos en el ADN determinan la información genética que codifica para las proteínas.

**Similitud Genética**: Los seres vivos comparten un alto grado de similitud genética en sus secuencias de ADN. Esta similitud se refleja en las secuencias de proteínas que producen. Las pequeñas diferencias genéticas entre especies pueden dar lugar a las diversas características y funciones biológicas que los hacen únicos, pero la mayor parte del ADN es común entre organismos cercanamente relacionados.
""")

# Incluir una imagen relacionada con genética
st.sidebar.image("https://images.my.labster.com/632b09d9-12b5-4bc1-937e-86b46c085d23/Codon_circle.es_ES.png")

# Entrada del usuario: Secuencia 1
st.subheader("Secuencia 1")
dna_input1 = st.text_area("Introduce la primera secuencia de ADN:", height=100)

# Entrada del usuario: Secuencia 2
st.subheader("Secuencia 2")
dna_input2 = st.text_area("Introduce la segunda secuencia de ADN:", height=100)

# Sidebar para seleccionar qué visualizar
view_option = st.sidebar.selectbox(
    "Selecciona qué visualizar:",
    ("Nucleótidos", "Proteínas")
)

# Procesamiento al pulsar el botón
if st.button("Procesar Secuencias"):
    results = []
    if dna_input1 or dna_input2:  # Verifica que al menos una secuencia sea ingresada
        for idx, dna_input in enumerate([dna_input1, dna_input2], start=1):
            if dna_input:
                # Quitar espacios y convertir a mayúsculas
                dna_sequence = dna_input.replace(" ", "").upper()

                # Validar que solo contenga caracteres válidos (A, T, C, G)
                if all(base in "ATCG" for base in dna_sequence):
                    protein_sequence = translate_dna_to_protein(dna_sequence)
                    protein_counts = count_proteins(protein_sequence)

                    # Conteo de nucleótidos
                    nucleotide_counts = count_nucleotides(dna_sequence)

                    # Traducción a proteínas
                    protein_sequence = translate_dna_to_protein(dna_sequence)
                    protein_counts = count_proteins(protein_sequence)

                    # Guarda los resultados en 'results'
                    results.append((nucleotide_counts, protein_counts))
                
                else:
                    st.error(f"La secuencia {idx} contiene caracteres inválidos. Por favor, introduce solo A, T, C y G.")
            else:
                st.warning(f"No se ingresó la Secuencia {idx}.")
    else:
        st.error("Por favor, introduce al menos una secuencia de ADN.")

        # Visualización según la opción seleccionada
        if results:
            for idx, (nucleotide_counts, protein_counts) in enumerate(results, start=1):
                if view_option == "Nucleótidos":
                    st.subheader(f"Nucleótidos en la Secuencia {idx}")
                    st.markdown("**Conteo de Nucleótidos**")
                    st.write(dict(nucleotide_counts))

                    st.markdown("**Gráfica de Nucleótidos**")
                    fig = plot_single_bar(
                        nucleotide_counts,
                        "Conteo de Nucleótidos",
                        "Nucleótidos",
                        "Cantidad",
                        color="lightblue"
                    )
                    st.pyplot(fig)

                elif view_option == "Proteínas":
                    st.subheader(f"Proteínas en la Secuencia {idx}")
                    st.markdown("**Conteo de Proteínas**")
                    st.write(dict(protein_counts))

                    st.markdown("**Gráfica de Proteínas**")
                    fig = plot_single_bar(
                        protein_counts,
                        "Conteo de Proteínas",
                        "Proteínas",
                        "Cantidad",
                        color="orange"
                    )
                    st.pyplot(fig)
        
         # Gráficas apiladas (si ambas secuencias son válidas)
        if len(results) == 2:
            nucleotides_labels = sorted(set(results[0][0].keys()).union(results[1][0].keys()))
            st.markdown("### Gráfico de Barras Apiladas - Nucleótidos")
            nucleotide_fig = plot_stacked_bar(
                results[0][0], results[1][0], nucleotides_labels, 
                "Comparación de Nucleótidos", "Nucleótidos", "Cantidad", 
                color1="lightgreen", color2="lightblue"
            )
            st.write("Figura de nucléotidos:", nucleotide_fig)
            st.pyplot(nucleotide_fig)

            proteins_labels = sorted(set(results[0][1].keys()).union(results[1][1].keys()))
            st.markdown("### Gráfico de Barras Apiladas - Proteínas")
            protein_fig = plot_stacked_bar(
                results[0][1], results[1][1], proteins_labels, 
                "Comparación de Proteínas", "Proteínas", "Cantidad", 
                color1="skyblue", color2="orange"
            )
            st.write("Figura de proteínas:", protein_fig)
            st.pyplot(protein_fig)
        else:
            st.warning("Se necesitan dos secuencias válidas para generar gráficos apilados.")
