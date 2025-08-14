# mview

A nimble metagenomic visualization tool for exploring genome assemblies.

## Description

Mview provides a customizable R-based shiny application for visualizing genomic data. The tool is designed for researchers working with metagenomic assemblies.

Mview is built around four key components that work together to create flexible genomic visualizations:

- **Configurations**: Define data sources and how to visualize them.

- **Views**: Define a set profiles to display. Each configuration defines one or more views. 

- **Profiles**: Individual visualization components like gene tracks, alignment displays, and coordinate axes.

- **Tabs**: Contigs, genomes, genes, read alignments and other tables.

Users visualize their data by creating configuration files that point to their tables while using existing profiles. The architecture supports the implementation of new profiles.

## Prerequisites

### Installation

Install required packages in R:

```r
install.packages(c("shiny", "DT", "plotly", "ggplot2", "shinyjs", "shinyjqui"))
```

### Installing alntools (required for the read alignment profile)

Alntools is a C++ program that allows for efficient querying of alignments for real-time visualization in mview.

Install **alntools** from GitHub:

```bash
# Clone and build from source
git clone https://github.com/eitanyaffe/alntools.git
cd alntools
make install
```

Add the local bin directory to your path if needed:
```
export PATH="$HOME/.local/bin:$PATH"
```

Install R package:
```r
install.packages("Rcpp")
```

## Quick Start

### 1. Start R and launch mview

```r
# Start R in the mview directory
source("mview.r")
```

### 2. Load a configuration

Once the app starts, you need to load a configuration. In the R console, type:

```r
# Load the minimal example configuration
rl("minimal")
```

The `rl()` function loads a configuration defined in `configs/minimal/minimal_cfg.r`

## Using mview

Application layout:

- **Central region**: Plot region with genomic profiles (also called tracks).
- **Left panel**: Information about the display window and cursor location.
- **Right panel**: Window where user-defined profiles parameters can be customized.
- **Bottom tabs**: Tables (states, contigs, genomes, genes, alignments).

#### Selecting the viewing region

1. **Select contigs**: Use State, Contigs or Genomes tab to select one or more contigs
2. **Zoom**: Use mouse to select regions in the visualization
3. **Lock zoom**: Press Shift+Z to set the selected region as zoom
4. **Navigation**: Use Shift+Backspace for undo, press Help button for all shortcuts

#### Tables

The bottom panel contains these tabs:
- **States**: Manage saved visualization states 
- **Contigs**: Table for selecting individual contigs
- **Genomes**: Table for selecting reconstructed genomes
- **Contig Map**: Mapping between contigs and genomes
- **Selected Contigs**: Displays currently selected contigs
- **Options**: Display settings

#### Read alignment tab

The Alignments tab shows all alignment of a selected read. A read can be selected by clicking on it in the alignment profile (if multiple alignment profiles are shown only the bottom profile can be clicked). The table shows all alignments for the selected read. Each alignment entry includes a "Go" button that navigates the view to the specific genomic location of that alignment, automatically setting the zoom and contig selection.

#### Gene tab

The Genes tab details all genes in view. Click on any gene in the table to select it for navigation. The "Navigate to Gene" button zooms the visualization to the selected gene's location with appropriate padding, making it easy to examine genes of interest in their genomic context.

#### Profile parameters
The collapsible right panel contains profile-specific settings that are applied in real-time to customize visualization appearance and behavior. The specific parameters shown here can be controlled through the configuration file.

### Saving and loading states

**Initialize state table**: Create a new table ("New Table").

**To save current state**:
1. Navigate to States tab
2. Enter description in text field  
3. Click "Add Current State"
4. Click "Save Table" to save the state table

**To restore a state**:
1. Load an existing state table by clicking "Load Table"
2. Select state from table
3. Click "Go to Selected State"

## Working with your data

See [data_setup.md](docs/data_setup.md) for detailed instructions on preparing data files and configuration.

## Customizing views and profiles

**Views**: Create custom views by copying existing view files and modifying which profiles to display and their parameters. Each view defines a specific layout of gene tracks, read alignments and other profiles.

**Profiles**: Building custom profiles is straightforward - use existing profile code from the `profiles/` directory as a template and modify the plotting functions to handle your specific data types and visualization needs.
