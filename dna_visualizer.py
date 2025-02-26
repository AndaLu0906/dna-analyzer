import tkinter as tk
from tkinter import ttk, scrolledtext
import random
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from DNASequenceAnalyzer import DNASequenceAnalyzer  # Assuming your original class is in a separate file

class DNAAnalyzerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("DNA Sequence Analyzer")
        self.root.geometry("800x600")
        self.analyzer = DNASequenceAnalyzer()

        # Create notebook (tabs)
        self.notebook = ttk.Notebook(root)
        self.notebook.pack(pady=10, fill="both", expand=True)

        # Tabs
        self.main_tab = ttk.Frame(self.notebook)
        self.results_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.main_tab, text="Input")
        self.notebook.add(self.results_tab, text="Results")

        # Main tab setup
        self.setup_main_tab()

        # Results tab setup (initially empty)
        self.setup_results_tab()

    def setup_main_tab(self):
        # Input frame
        input_frame = ttk.LabelFrame(self.main_tab, text="DNA Sequence Input")
        input_frame.pack(pady=10, padx=10, fill="x")

        ttk.Label(input_frame, text="Enter DNA Sequence (or leave blank for random):").pack(pady=5)
        self.sequence_entry = scrolledtext.ScrolledText(input_frame, height=5, width=50)
        self.sequence_entry.pack(pady=5)

        ttk.Label(input_frame, text="Sequence Length (for random generation):").pack(pady=5)
        self.length_entry = ttk.Entry(input_frame)
        self.length_entry.insert(0, "2000")
        self.length_entry.pack(pady=5)

        # Analyze button
        ttk.Button(input_frame, text="Analyze", command=self.run_analysis).pack(pady=10)

    def setup_results_tab(self):
        self.results_text = scrolledtext.ScrolledText(self.results_tab, height=20, width=80)
        self.results_text.pack(pady=10, padx=10, fill="both", expand=True)

    def run_analysis(self):
        # Clear previous results
        self.results_text.delete(1.0, tk.END)

        # Get sequence
        sequence = self.sequence_entry.get(1.0, tk.END).strip()
        if not sequence:
            try:
                length = int(self.length_entry.get())
                sequence = self.analyzer.generate_random_dna(length)
            except ValueError:
                self.results_text.insert(tk.END, "Invalid length entered!")
                return

        if not self.analyzer.validate_dna_sequence(sequence):
            self.results_text.insert(tk.END, "Invalid DNA sequence! Use only A, T, G, C.")
            return

        # Basic analysis
        self.results_text.insert(tk.END, f"Sequence Length: {len(sequence)}\n")
        self.results_text.insert(tk.END, f"First 60 nucleotides: {sequence[:60]}...\n\n")

        # Nucleotide frequency (with matplotlib bar chart)
        freq = self.analyzer.nucleotide_frequency(sequence)
        self.show_bar_chart(freq, "Nucleotide Frequency", "Nucleotide", "Count")

        # GC content
        gc = self.analyzer.gc_content(sequence)
        self.results_text.insert(tk.END, f"GC Content: {gc:.2f}%\n\n")

        # GC content variation
        positions, gc_values = self.analyzer.gc_content_regions(sequence)
        self.show_line_plot(positions, gc_values, "GC Content Variation", "Position", "GC %")

        # Find ORFs (showing first one as example)
        orfs = self.analyzer.find_open_reading_frames(sequence)
        if orfs:
            strand, frame, protein = orfs[0]
            self.results_text.insert(tk.END, f"Example ORF ({strand} strand, frame {frame}):\n")
            self.results_text.insert(tk.END, f"Protein: {protein[:30]}... (length: {len(protein)})\n")
            analysis = self.analyzer.analyze_protein(protein)
            self.results_text.insert(tk.END, f"Hydrophobic Ratio: {analysis['hydrophobic_ratio']:.2f}\n")
            self.results_text.insert(tk.END, f"Net Charge: {analysis['charge']}\n\n")

    def show_bar_chart(self, data, title, xlabel, ylabel):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.bar(data.keys(), data.values())
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        self.display_plot(fig)

    def show_line_plot(self, x, y, title, xlabel, ylabel):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(x, y, '-o', markersize=4)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        self.display_plot(fig)

    def display_plot(self, fig):
        new_window = tk.Toplevel(self.root)
        new_window.title("Chart")
        canvas = FigureCanvasTkAgg(fig, master=new_window)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

if __name__ == "__main__":
    root = tk.Tk()
    app = DNAAnalyzerGUI(root)
    root.mainloop()
