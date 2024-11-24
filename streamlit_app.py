import streamlit as st
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

# Configuración de Streamlit
st.title("Dashboard de Traducción de ADN a Proteínas")
st.markdown("""
Este dashboard permite introducir una secuencia de nucleótidos (ADN) para traducirla 
en proteínas y contar la cantidad de cada una presente en la secuencia.
""")

# Entrada del usuario
dna_input = st.text_area("Introduce la secuencia de ADN:", height=200)

# Procesamiento al pulsar el botón
if st.button("Procesar"):
    if dna_input:
        # Quitar espacios y convertir a mayúsculas
        dna_sequence = dna_input.replace(" ", "").upper()

        # Validar que solo contenga caracteres válidos (A, T, C, G)
        if all(base in "ATCG" for base in dna_sequence):
            protein_sequence = translate_dna_to_protein(dna_sequence)
            protein_counts = count_proteins(protein_sequence)

            # Mostrar resultados
            st.subheader("Secuencia de proteínas traducida")
            st.text(protein_sequence)

            st.subheader("Conteo de proteínas")
            st.write(dict(protein_counts))
        else:
            st.error("La secuencia contiene caracteres inválidos. Por favor, introduce solo A, T, C y G.")
    else:
        st.error("Por favor, introduce una secuencia de ADN.")
