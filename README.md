# DNA Analyzer

A Python-based tool for analyzing DNA sequences with a graphical user interface (GUI) built using Tkinter and Matplotlib. This project provides functionalities to generate, validate, and analyze DNA sequences, including nucleotide frequency, GC content, open reading frames (ORFs), and protein properties.

## Features
- **Random DNA Generation**: Generate random DNA sequences of specified lengths.
- **Sequence Validation**: Ensures input sequences contain only A, T, G, C.
- **Nucleotide Frequency**: Visualize the distribution of A, T, G, and C with bar charts.
- **GC Content Analysis**: Calculate and plot GC content across sequence regions.
- **ORF Detection**: Identify potential open reading frames in forward and reverse strands.
- **Protein Analysis**: Analyze translated protein sequences for hydrophobicity, charge, and amino acid composition.
- **GUI Interface**: User-friendly input and results display using Tkinter, with plots in separate windows via Matplotlib.

## Files
- `dna_visualizer.py`: The main script with the Tkinter GUI.
- `DNASequenceAnalyzer.py`: The core class handling DNA sequence analysis logic.

## Requirements
- Python 3.13+
- Tkinter (included with Python, requires Tcl/Tk)
- Matplotlib (`pip install matplotlib`)

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/AndaLu0906/dna_analyzer.git
   cd dna_analyzer≈
