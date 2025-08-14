# mview

A nimble metagenomic visualization tool for exploring genome assemblies.

## Description

Mview provides a customizable R-based shiny application for visualizing genomic data. The tool is designed for researchers working with reference-based and reconstructed metagenomic assemblies.

Mview is built around four key components that work together to create flexible genomic visualizations:

- **Configurations**: Define data sources, assemblies, available tabs, and which views to load. These files specify how to read your genomic data tables and register visualization components.

- **Views**: Define a set profiles to display. Configurations define one or more views. 

- **Profiles**: Individual visualization components like gene tracks, alignment displays, and coordinate axes. Each profile handles a specific type of genomic data visualization.

- **Tabs**: Optional bottom panel components that provide detailed data tables and interactive features for genes, alignments, and other data types.

Users visualize their data by creating configuration files that point to their tables while using existing profiles. Advanced users can implement custom profiles to create new visualizations.

## Prerequisites

### Installation

Install all required packages at once:

```r
install.packages(c("shiny", "DT", "plotly", "ggplot2", "shinyjs", "shinyjqui"))
```

### Installing alntools (Optional)

Alntools is a C++ program that allows for efficient querying of alignments for real-time visualization.

To visualize read alignments, you need to install **alntools** from GitHub:

```bash
# Clone and build from source
git clone https://github.com/eitanyaffe/alntools.git
cd alntools
make
# Copy binary to user's local bin directory
mkdir -p ~/bin
cp bin/alntools ~/bin/
export PATH="$HOME/bin:$PATH"
```

Install required R package:
```r
install.packages(c("Rcpp"))
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

### 3. Basic usage of mview

Application layout:

- **Left**: Information about display region.
- **Middle**: Profiles plots.
- **Right**: Profile parameter window.
- **Bottom**: Tables in tabs.

#### Selecting the viewing region

1. **Select contigs**: Use Contigs or Genomes tab to select one or more contigs
2. **Zoom**: Use mouse to select regions in the visualization
3. **Lock zoom**: Press Shift+Z to set the selected region as zoom
4. **Navigation**: Use Shift+Backspace for undo, press Help button for all shortcuts

**State management**: Use the States tab to save current view (contigs + zoom) and restore saved states.

#### Tables

The bottom panel contains these tabs:
- **States**: Manage saved visualization states 
- **Contigs**: Table for selecting individual contigs
- **Genomes**: Table for selecting reconstructed genomes
- **Contig Map**: Mapping between contigs and genomes
- **Selected Contigs**: Displays currently selected contigs
- **Options**: Display settings

**Additinal tabs** (optional; shown only if registered in configuration):
- **Genes**: Gene table
- **Alignments**: Read alignment table

#### Profile parameters
The collapsible right panel contains profile-specific settings that are applied in real-time to customize visualization appearance and behavior.


## Saving and loading state

**Initialize state table**: Create a new table ("New Table") or load an existing state table ("Load Table").

**To save current state**:
1. Navigate to States tab
2. Enter description in text field  
3. Click "Add Current State"

**To restore a state**:
1. Select state from table
2. Click "Go to Selected State"

**Additional actions**: "Delete Selected State", "Save Table" (to file), "Load Table" (from file)

## Working with your data

See [data_setup.md](docs/data_setup.md) for detailed instructions on preparing data files and configuration.

## Customizing views and profiles

**Views**: Create custom views by copying existing view files and modifying which profiles to display and their parameters. Each view defines a specific layout of gene tracks, alignment displays, and visualization settings.

**Profiles**: Building custom profiles is straightforward - copy existing profile code from the `profiles/` directory and modify the plotting functions to handle your specific data types and visualization needs.
