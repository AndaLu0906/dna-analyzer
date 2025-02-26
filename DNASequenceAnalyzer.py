import random
from collections import Counter

class DNASequenceAnalyzer:
    """A tool for analyzing DNA sequences and simulating protein translation."""
    
    def __init__(self):
        # Genetic code (codon to amino acid mapping)
        self.genetic_code = {
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
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
            'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
        }
        
        # Properties of amino acids
        self.amino_acid_properties = {
            'A': {'hydrophobic': True, 'charge': 0, 'size': 'small'},
            'C': {'hydrophobic': True, 'charge': 0, 'size': 'small'},
            'D': {'hydrophobic': False, 'charge': -1, 'size': 'small'},
            'E': {'hydrophobic': False, 'charge': -1, 'size': 'medium'},
            'F': {'hydrophobic': True, 'charge': 0, 'size': 'large'},
            'G': {'hydrophobic': True, 'charge': 0, 'size': 'small'},
            'H': {'hydrophobic': False, 'charge': 1, 'size': 'medium'},
            'I': {'hydrophobic': True, 'charge': 0, 'size': 'large'},
            'K': {'hydrophobic': False, 'charge': 1, 'size': 'large'},
            'L': {'hydrophobic': True, 'charge': 0, 'size': 'large'},
            'M': {'hydrophobic': True, 'charge': 0, 'size': 'large'},
            'N': {'hydrophobic': False, 'charge': 0, 'size': 'small'},
            'P': {'hydrophobic': False, 'charge': 0, 'size': 'small'},
            'Q': {'hydrophobic': False, 'charge': 0, 'size': 'medium'},
            'R': {'hydrophobic': False, 'charge': 1, 'size': 'large'},
            'S': {'hydrophobic': False, 'charge': 0, 'size': 'small'},
            'T': {'hydrophobic': False, 'charge': 0, 'size': 'small'},
            'V': {'hydrophobic': True, 'charge': 0, 'size': 'medium'},
            'W': {'hydrophobic': True, 'charge': 0, 'size': 'large'},
            'Y': {'hydrophobic': True, 'charge': 0, 'size': 'large'},
            '*': {'hydrophobic': False, 'charge': 0, 'size': 'none'},
        }
        
    def validate_dna_sequence(self, sequence):
        """Validate that a sequence contains only A, T, G, C."""
        valid_nucleotides = set('ATGC')
        return all(nucleotide in valid_nucleotides for nucleotide in sequence.upper())
    
    def generate_random_dna(self, length=1000):
        """Generate a random DNA sequence of specified length."""
        return ''.join(random.choice('ATGC') for _ in range(length))
    
    def nucleotide_frequency(self, sequence):
        """Calculate the frequency of each nucleotide in the sequence."""
        sequence = sequence.upper()
        counts = Counter(sequence)
        return {base: counts.get(base, 0) for base in 'ATGC'}
    
    def gc_content(self, sequence):
        """Calculate the GC content percentage."""
        sequence = sequence.upper()
        g_count = sequence.count('G')
        c_count = sequence.count('C')
        return (g_count + c_count) / len(sequence) * 100 if len(sequence) > 0 else 0
    
    def find_motifs(self, sequence, motif):
        """Find all occurrences of a specific motif in the sequence."""
        sequence = sequence.upper()
        motif = motif.upper()
        positions = []
        for i in range(len(sequence) - len(motif) + 1):
            if sequence[i:i+len(motif)] == motif:
                positions.append(i)
        return positions
    
    def translate_dna(self, sequence, start_pos=0):
        """
        Translate a DNA sequence into a protein sequence.
        Starts from the given position and goes until a stop codon or end of sequence.
        """
        sequence = sequence.upper()
        protein = ""
        
        # Find the first start codon (ATG)
        start_index = sequence.find('ATG', start_pos)
        if start_index == -1:  # No start codon found
            return ""
        
        # Translate codons to amino acids
        for i in range(start_index, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if len(codon) < 3:  # Incomplete codon at the end
                break
                
            amino_acid = self.genetic_code.get(codon, 'X')  # 'X' for unknown codons
            if amino_acid == '*':  # Stop codon
                break
                
            protein += amino_acid
            
        return protein
    
    def find_open_reading_frames(self, sequence, min_length=30):
        """Find all potential open reading frames in the sequence."""
        sequence = sequence.upper()
        orfs = []
        
        # Check all six reading frames (3 forward, 3 reverse complement)
        for frame in range(3):
            # Forward sequence
            protein = self.translate_dna(sequence, frame)
            if len(protein) >= min_length / 3:  # Convert min_length from nucleotides to amino acids
                orfs.append(("forward", frame, protein))
            
            # Reverse complement sequence
            rev_comp = self.reverse_complement(sequence)
            protein = self.translate_dna(rev_comp, frame)
            if len(protein) >= min_length / 3:
                orfs.append(("reverse", frame, protein))
        
        return orfs
    
    def reverse_complement(self, sequence):
        """Get the reverse complement of a DNA sequence."""
        complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(complement_map.get(base, 'N') for base in reversed(sequence.upper()))
    
    def analyze_protein(self, protein_sequence):
        """Analyze properties of a protein sequence."""
        if not protein_sequence:
            return {
                "length": 0,
                "hydrophobic_ratio": 0,
                "charge": 0,
                "amino_acid_counts": {}
            }
            
        # Count amino acids
        aa_counts = Counter(protein_sequence)
        
        # Calculate hydrophobic ratio
        hydrophobic_count = sum(1 for aa in protein_sequence if self.amino_acid_properties.get(aa, {}).get('hydrophobic', False))
        hydrophobic_ratio = hydrophobic_count / len(protein_sequence) if len(protein_sequence) > 0 else 0
        
        # Calculate net charge
        total_charge = sum(self.amino_acid_properties.get(aa, {}).get('charge', 0) for aa in protein_sequence)
        
        return {
            "length": len(protein_sequence),
            "hydrophobic_ratio": hydrophobic_ratio,
            "charge": total_charge,
            "amino_acid_counts": dict(aa_counts)
        }
    
    def text_histogram(self, data, width=50, label="Value"):
        """Generate a text-based histogram for displaying data."""
        if not data:
            return "No data to display"
            
        result = []
        max_value = max(data.values())
        max_label_len = max(len(str(k)) for k in data.keys())
        
        for key, value in sorted(data.items()):
            bar_length = int((value / max_value) * width) if max_value > 0 else 0
            bar = "█" * bar_length
            result.append(f"{key:{max_label_len}} | {bar} {value}")
            
        return "\n".join(result)
    
    def gc_content_regions(self, sequence, window_size=100, step_size=50):
        """
        Calculate GC content across regions of the sequence.
        Returns positions and corresponding GC percentages.
        """
        positions = []
        gc_values = []
        
        for i in range(0, len(sequence) - window_size + 1, step_size):
            window = sequence[i:i+window_size]
            gc = self.gc_content(window)
            positions.append(i + window_size // 2)
            gc_values.append(gc)
            
        return positions, gc_values
    
    def text_plot(self, x_values, y_values, width=80, height=20, title="", x_label="", y_label=""):
        """Generate a simple text-based plot."""
        if not x_values or not y_values or len(x_values) != len(y_values):
            return "Invalid data for plotting"
            
        # Calculate min and max for scaling
        min_x, max_x = min(x_values), max(x_values)
        min_y, max_y = min(y_values), max(y_values)
        
        # Handle edge cases
        if min_x == max_x:
            max_x = min_x + 1
        if min_y == max_y:
            max_y = min_y + 1
            
        # Create the empty plot grid
        grid = [[' ' for _ in range(width)] for _ in range(height)]
        
        # Plot the points
        for x, y in zip(x_values, y_values):
            # Scale to grid coordinates
            grid_x = int((x - min_x) / (max_x - min_x) * (width - 1))
            grid_y = height - 1 - int((y - min_y) / (max_y - min_y) * (height - 1))
            
            # Ensure within bounds
            grid_x = max(0, min(grid_x, width - 1))
            grid_y = max(0, min(grid_y, height - 1))
            
            grid[grid_y][grid_x] = '●'
            
        # Convert grid to string
        result = [title.center(width)]
        for row in grid:
            result.append(''.join(row))
        result.append(x_label.center(width))
        
        # Add y-axis labels
        temp_result = []
        if y_label:
            y_labels = [f"{min_y:.1f}", f"{(min_y + max_y) / 2:.1f}", f"{max_y:.1f}"]
            for i, line in enumerate(result):
                if i == 0:
                    temp_result.append(f"{' ' * (len(y_labels[2]) + 2)}{line}")
                elif i == 1:
                    temp_result.append(f"{y_labels[2]} {line}")
                elif i == height // 2:
                    temp_result.append(f"{y_labels[1]} {line}")
                elif i == height:
                    temp_result.append(f"{y_labels[0]} {line}")
                else:
                    temp_result.append(f"{' ' * (len(y_labels[2]) + 2)}{line}")
            result = temp_result
            
        return '\n'.join(result)
    
    def calculate_hydrophobicity_profile(self, protein_sequence, window_size=7):
        """Calculate the hydrophobicity profile of a protein sequence."""
        if not protein_sequence or len(protein_sequence) < window_size:
            return [], []
            
        # Hydrophobicity values (Kyte-Doolittle scale)
        hydrophobicity = {
            'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
            'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
            'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
            'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
        }
        
        # Calculate average hydrophobicity in sliding windows
        hydrophobicity_values = []
        positions = []
        
        for i in range(len(protein_sequence) - window_size + 1):
            window = protein_sequence[i:i+window_size]
            avg_hydrophobicity = sum(hydrophobicity.get(aa, 0) for aa in window) / window_size
            hydrophobicity_values.append(avg_hydrophobicity)
            positions.append(i + window_size // 2)
        
        return positions, hydrophobicity_values
    
    def identify_domains(self, protein_sequence):
        """
        Identify potential protein domains based on hydrophobicity patterns.
        This is a simplified approach for demonstration purposes.
        """
        if not protein_sequence or len(protein_sequence) < 15:
            return []
            
        # Calculate hydrophobicity profile
        positions, hydrophobicity = self.calculate_hydrophobicity_profile(protein_sequence, window_size=9)
        
        # Identify regions with consistent hydrophobicity patterns
        domains = []
        in_domain = False
        start_pos = 0
        domain_type = ""
        
        for i, h_value in enumerate(hydrophobicity):
            # Hydrophobic domain
            if h_value > 1.0 and not in_domain:
                in_domain = True
                start_pos = positions[i]
                domain_type = "hydrophobic"
            # Hydrophilic domain
            elif h_value < -1.0 and not in_domain:
                in_domain = True
                start_pos = positions[i]
                domain_type = "hydrophilic"
            # End of domain
            elif in_domain and ((domain_type == "hydrophobic" and h_value < 0) or
                              (domain_type == "hydrophilic" and h_value > 0)):
                end_pos = positions[i]
                if end_pos - start_pos >= 9:  # Minimum domain length
                    domains.append({
                        "type": domain_type,
                        "start": start_pos,
                        "end": end_pos,
                        "length": end_pos - start_pos + 1
                    })
                in_domain = False
        
        # Check for domain at the end of the sequence
        if in_domain and positions:
            end_pos = positions[-1]
            if end_pos - start_pos >= 9:
                domains.append({
                    "type": domain_type,
                    "start": start_pos,
                    "end": end_pos,
                    "length": end_pos - start_pos + 1
                })
        
        return domains
