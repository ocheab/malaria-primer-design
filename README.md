# malariaPrimerDesigner

**malariaPrimerDesigner** is an R package and Shiny web application developed by the *Centre for Malaria and Other Tropical Diseases Care (CEMTROD)*, University of Ilorin Teaching Hospital, Nigeria. The tool provides an intuitive interface for designing high-quality PCR primers for diagnostic markers across multiple *Plasmodium* species.

---

## 🔬 Overview

Malaria diagnosis often relies on molecular markers such as **MSP1**, **HRP2**, and **18S rRNA**. Designing efficient and specific primers for these markers can be time-consuming and error-prone. `malariaPrimerDesigner` simplifies this process by offering a robust graphical interface and automated filtering system based on best practices in primer design.

The tool enables:

- 📂 **Built-in marker support**: Design primers for validated malaria diagnostic genes using preloaded species-specific reference sequences.
- 📁 **Custom sequence uploads**: Upload and design primers from your own FASTA sequences.
- 🧬 **Flexible filtering options**:
  - Adjustable primer **length** range.
  - User-defined **GC content** and **melting temperature (Tm)** range.
  - Detection of **GC clamp** at 3’ ends for stable annealing.
  - Filtering of primers with potential **hairpin loops** or **cross-dimers**.
- 📊 **Visual diagnostics**: Histogram plots of Tm and GC content distributions help assess primer quality.
- 📤 **Export functionality**:
  - Download the full table of filtered primers as a `.csv` file.
  - Select and export primer pairs in **FASTA** format.

---

## 🧪 Features at a Glance

| Feature                  | Description                                                 |
|--------------------------|-------------------------------------------------------------|
| Primer Design            | For MSP1, HRP2, 18S rRNA and user-uploaded sequences        |
| GC Clamp Detection       | Ensures strong binding at 3' ends                           |
| Hairpin & Dimer Checks   | Screens for structural instability                          |
| Custom Filtering         | Tm, GC content, and primer length thresholds                |
| Multiple Species Support | *P. falciparum*, *P. vivax*, *P. malariae*, etc.            |
| FASTA & CSV Export       | Output selected primers for downstream lab use              |

---

## 📦 Installation

The package is currently hosted on GitHub.

To install:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install malariaPrimerDesigner from GitHub
devtools::install_github("ocheab/malaria-primer-design", INSTALL_opts = "--no-staged-install")
