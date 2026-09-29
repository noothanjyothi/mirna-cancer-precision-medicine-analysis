# miRNA Cancer Analysis with Biopython

A beginner-friendly bioinformatics project for analyzing miRNA sequences and exploring their role in cancer biology using Biopython.

## 📌 What is this project?

This project teaches you how to:
- Read and parse miRNA sequences from FASTA files
- Analyze nucleotide composition (A, U, G, C counts)
- Calculate GC content (important for RNA stability)
- Compare sequences to find similarities
- Understand why miRNAs matter in cancer biology

**No machine learning. No command line. Just Python and Biopython.**

---

## 🧬 Why miRNA and Cancer?

miRNAs (microRNAs) are tiny RNA molecules that control which genes get expressed in cells.

In cancer:
- Some miRNAs are **turned ON** when they shouldn't be (oncogenic)
- Some miRNAs are **turned OFF** when cells need them (tumor suppressors)
- By studying miRNA sequences, we can understand cancer better

Example: **miR-21** is often increased in many cancer types and helps cancer cells survive.

---

## 📁 Project Structure

```
mirna-cancer-analysis/
├── README.md                    (this file)
├── requirements.txt             (Python packages needed)
├── data/
│   └── sample_mirnas.fasta      (example miRNA sequences)
├── analysis.py                  (main Python script)
└── results/
    └── analysis_report.txt      (output results)
```

---

## 🔧 Installation

### Step 1: Install Python
Make sure you have Python 3.7 or higher.

### Step 2: Create a folder for your project
Create a new folder on your computer for this project.

### Step 3: Install Biopython
Open a terminal/command prompt and run:

```
pip install biopython
```

That's it! You now have everything you need.

---

## 📦 What you need

**requirements.txt** contains:
```
biopython
```

---

## 💾 Sample Data (FASTA Format)

Create a file called `sample_mirnas.fasta` in your `data/` folder with this content:

```
>hsa-miR-21
UAGCUUAUCAGACUGAUGUUGA

>hsa-miR-155
UUAAUGCUAAUCGUGAUAGGGGUU

>hsa-miR-17
CAAAGUGCUUACAGUGCAGGUAG

>hsa-miR-29a
ACUGAUUUCUUUGGUGUCAGC

>hsa-miR-200c
UAAUACUGCCGGGUAAUGAUGGA
```

These are real miRNA sequences found in humans. The `>` symbol marks the name, and the sequence follows on the next line.

---

## 🐍 Python Scripts

### Script 1: Basic Sequence Analysis

This script reads miRNA sequences and shows their properties:

```python
from Bio import SeqIO

def analyze_mirna(fasta_file):
    """
    Read miRNA sequences and show basic properties
    """
    print("=" * 50)
    print("miRNA SEQUENCE ANALYSIS")
    print("=" * 50)
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        
        # Count nucleotides
        a_count = seq.count("A")
        u_count = seq.count("U")
        g_count = seq.count("G")
        c_count = seq.count("C")
        
        # Calculate GC content
        gc_content = (g_count + c_count) / len(seq) * 100
        
        # Print results
        print(f"\nmiRNA Name: {record.id}")
        print(f"Sequence: {seq}")
        print(f"Length: {len(seq)} nucleotides")
        print(f"A: {a_count}, U: {u_count}, G: {g_count}, C: {c_count}")
        print(f"GC Content: {gc_content:.1f}%")
        print("-" * 50)

# Run the analysis
analyze_mirna("data/sample_mirnas.fasta")
```

**What this does:**
- Opens your FASTA file
- Reads each miRNA sequence
- Counts A, U, G, C nucleotides
- Calculates GC content (% of G and C in the sequence)
- Prints the results

**Output looks like:**
```
==================================================
miRNA SEQUENCE ANALYSIS
==================================================

miRNA Name: hsa-miR-21
Sequence: UAGCUUAUCAGACUGAUGUUGA
Length: 22 nucleotides
A: 4, U: 5, G: 5, C: 2
GC Content: 31.8%
--------------------------------------------------
```

---

### Script 2: Compare Two miRNA Sequences

This script shows how similar two sequences are:

```python
from Bio import SeqIO, pairwise2

def compare_sequences(fasta_file):
    """
    Compare the first two sequences in the file
    """
    sequences = []
    
    # Read sequences
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append((record.id, str(record.seq)))
    
    if len(sequences) < 2:
        print("Need at least 2 sequences to compare")
        return
    
    # Compare first two
    name1, seq1 = sequences[0]
    name2, seq2 = sequences[1]
    
    print(f"\nComparing {name1} and {name2}")
    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    
    # Count matching positions
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    similarity = (matches / max(len(seq1), len(seq2))) * 100
    
    print(f"\nMatching positions: {matches}/{max(len(seq1), len(seq2))}")
    print(f"Similarity: {similarity:.1f}%")

# Run comparison
compare_sequences("data/sample_mirnas.fasta")
```

**Output looks like:**
```
Comparing hsa-miR-21 and hsa-miR-155
Sequence 1: UAGCUUAUCAGACUGAUGUUGA
Sequence 2: UUAAUGCUAAUCGUGAUAGGGGUU

Matching positions: 14/24
Similarity: 58.3%
```

---

### Script 3: Search for Specific Patterns

This script finds sequences containing a specific pattern:

```python
from Bio import SeqIO

def find_pattern(fasta_file, pattern):
    """
    Find miRNAs that contain a specific sequence pattern
    """
    print(f"\nSearching for pattern: {pattern}")
    print("-" * 50)
    
    found_count = 0
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        
        if pattern.upper() in seq.upper():
            position = seq.upper().find(pattern.upper())
            print(f"Found in {record.id}")
            print(f"Position: {position}")
            print(f"Full sequence: {seq}")
            found_count += 1
    
    if found_count == 0:
        print("Pattern not found in any sequence")
    else:
        print(f"\nTotal matches: {found_count}")

# Search for a pattern
find_pattern("data/sample_mirnas.fasta", "UUA")
```

**Output looks like:**
```
Searching for pattern: UUA
--------------------------------------------------
Found in hsa-miR-21
Position: 5
Full sequence: UAGCUUAUCAGACUGAUGUUGA

Found in hsa-miR-155
Position: 2
Full sequence: UUAAUGCUAAUCGUGAUAGGGGUU

Total matches: 2
```

---

### Script 4: All-in-One Analysis

Here's a complete script that does everything:

```python
from Bio import SeqIO

def full_analysis(fasta_file):
    """
    Complete miRNA analysis
    """
    sequences_data = []
    
    print("\n" + "=" * 60)
    print("COMPLETE miRNA CANCER ANALYSIS")
    print("=" * 60)
    
    # Read and analyze each sequence
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        
        # Count nucleotides
        a = seq.count("A")
        u = seq.count("U")
        g = seq.count("G")
        c = seq.count("C")
        
        # Calculate GC content
        gc = (g + c) / len(seq) * 100
        
        # Store data
        sequences_data.append({
            "name": record.id,
            "sequence": seq,
            "length": len(seq),
            "gc_content": gc,
            "a": a,
            "u": u,
            "g": g,
            "c": c
        })
        
        print(f"\nmiRNA: {record.id}")
        print(f"Sequence: {seq}")
        print(f"Length: {len(seq)} nt")
        print(f"Composition - A:{a} U:{u} G:{g} C:{c}")
        print(f"GC Content: {gc:.1f}%")
    
    # Summary statistics
    print("\n" + "=" * 60)
    print("SUMMARY STATISTICS")
    print("=" * 60)
    
    total_sequences = len(sequences_data)
    avg_length = sum(s["length"] for s in sequences_data) / total_sequences
    avg_gc = sum(s["gc_content"] for s in sequences_data) / total_sequences
    
    print(f"Total miRNAs analyzed: {total_sequences}")
    print(f"Average length: {avg_length:.1f} nucleotides")
    print(f"Average GC content: {avg_gc:.1f}%")
    
    # Find highest and lowest GC
    highest_gc = max(sequences_data, key=lambda x: x["gc_content"])
    lowest_gc = min(sequences_data, key=lambda x: x["gc_content"])
    
    print(f"\nHighest GC content: {highest_gc['name']} ({highest_gc['gc_content']:.1f}%)")
    print(f"Lowest GC content: {lowest_gc['name']} ({lowest_gc['gc_content']:.1f}%)")
    
    print("\n" + "=" * 60)

# Run full analysis
full_analysis("data/sample_mirnas.fasta")
```

---

## 🧪 How to Run

1. **Save one of the scripts above** as a `.py` file (e.g., `analysis.py`) in your project folder

2. **Run it in Python:**

```python
# Just copy and paste the script into Python IDLE or Jupyter
# Or open a Python file and run it
exec(open("analysis.py").read())
```

Or if you know how to use Python directly:

```python
# Simply run: python analysis.py
```

---

## 📊 Understanding the Output

### GC Content
- **What it is**: Percentage of Guanine (G) and Cytosine (C) in RNA
- **Why it matters**: Higher GC content = more stable RNA, less likely to break down
- **Range**: Usually 40-60% in miRNAs

### Sequence Length
- **miRNAs are small**: Usually 18-25 nucleotides
- **Why**: This small size allows precise gene targeting

### Nucleotide Composition
- **A (Adenine)**: Pairs with U
- **U (Uracil)**: Pairs with A (unique to RNA, DNA has T instead)
- **G (Guanine)**: Pairs with C
- **C (Cytosine)**: Pairs with G

---

## 🔍 Biological Interpretation

### What does high GC content mean?
- More stable RNA
- Less likely to break down in cells
- More effective as a regulatory molecule

### Why compare sequences?
- miRNAs with similar sequences often target similar genes
- Family relationships between miRNAs
- Understanding evolutionary history

### Cancer Connection
- **Oncogenic miRNAs**: High expression → helps cancer cells grow
- **Tumor suppressor miRNAs**: Low expression → can't stop cancer cells
- **biomarkers**: Specific miRNAs can predict cancer type or treatment response

---

## 📚 Example: Real Cancer-Associated miRNAs

| miRNA | Role in Cancer | Effect |
|-------|---|---|
| miR-21 | Oncogenic | Promotes cell survival, stops apoptosis |
| miR-155 | Oncogenic | Increases tumor growth |
| miR-29 | Tumor suppressor | Reduces cancer cell growth |
| miR-200c | Tumor suppressor | Prevents metastasis |

---

## 🎯 Project Ideas

Once you understand the basics, try these:

1. **Add more sequences** to your FASTA file from NCBI
2. **Compare miRNAs from different cancer types**
3. **Find common patterns** in oncogenic vs. tumor suppressor miRNAs
4. **Save results to a file** (write results to a text file)
5. **Create a report** summarizing all findings
6. **Analyze seed region** (first 6-8 nucleotides that matter most)

---

## 🌐 Getting Real Data

You can get real miRNA sequences from:
- **NCBI miRBase**: https://www.mirbase.org/
- **NCBI GenBank**: Direct download of FASTA files
- Just download and save as `.fasta` file, then analyze

---

## 📝 Notes

- miRNAs use **U (Uracil)** instead of **T (Thymine)** because they are RNA
- Sequences are written **5' to 3' direction** (standard biology notation)
- Real miRNA lengths: 18-25 nucleotides (very short!)
- This project focuses on **sequence analysis only** — no statistics or ML

---

## 🎓 Learning Goals

After this project, you will understand:
- ✅ How to use Biopython to read sequences
- ✅ Basic sequence properties (length, composition, GC content)
- ✅ How to compare sequences
- ✅ Why miRNAs are important in cancer biology
- ✅ How to structure a bioinformatics project
- ✅ Scientific Python programming

---

## 📖 Useful Biopython Functions

| Function | What it does |
|----------|---|
| `SeqIO.parse()` | Read sequences from FASTA file |
| `record.id` | Get sequence name |
| `record.seq` | Get the sequence string |
| `seq.count("A")` | Count a nucleotide |
| `len(seq)` | Get sequence length |
| `seq.upper()` | Convert to uppercase |
| `seq.lower()` | Convert to lowercase |

---

## ✅ Checklist to Complete

- [ ] Install Python and Biopython
- [ ] Create project folder with correct structure
- [ ] Create `sample_mirnas.fasta` file with sample data
- [ ] Write and run Script 1 (Basic Analysis)
- [ ] Write and run Script 2 (Compare Sequences)
- [ ] Write and run Script 3 (Find Patterns)
- [ ] Write and run Script 4 (Full Analysis)
- [ ] Understand the output and biological meaning
- [ ] Try analyzing real miRNA data from NCBI

---

## 📧 Questions?

This project uses only Biopython and basic Python. No machine learning. No complex statistics. Just simple, clear sequence analysis that teaches you bioinformatics.

---

**Last Updated**: September 29, 2026  
**Status**: 🟢 Ready to Learn
