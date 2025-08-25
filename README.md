# mview

A nimble metagenomic visualization tool for exploring genome assemblies.

## Description

Mview provides a customizable R-based shiny application for visualizing genomic data. The tool is designed for researchers working with metagenomic assemblies and features advanced mutation analysis including segregating sites detection and genomic distance visualization.

Mview is built around four key components that work together to create flexible genomic visualizations:

- **Configurations**: Define data sources and how to visualize them.

- **Views**: Define a set profiles to display. Each configuration defines one or more views. 

- **Profiles**: Individual visualization components like gene tracks, alignment displays, and coordinate axes.

- **Tabs**: Contigs, genomes, genes, read alignments and other tables.

Users visualize their data by creating configuration files that point to their tables while using existing profiles. The architecture supports the implementation of new profiles.

## Installation

### Installation Overview

mview is part of the makeshift toolkit for genomics analysis. The recommended setup organizes related tools in a shared directory structure:

```
makeshift/
└── tools/
    ├── mview/      # This visualization tool
    └── alntools/   # Optional alignment processing (for read alignment profiles)
```

This organization allows mview to seamlessly integrate with alntools for advanced alignment visualization features.

### Installing mview

Clone mview from GitHub:

```bash
# Create makeshift directory and clone mview
mkdir -p /path/to/your/makeshift/tools
cd /path/to/your/makeshift/tools
git clone https://github.com/eitanyaffe/mview.git
```

Install required packages in R:

```r
install.packages(c("shiny", "DT", "plotly", "ggplot2", "shinyjs", "shinyjqui", "seqinr"))
```

### Installing alntools (optional)

Install alntools if you wish to visualize read alignments. For optimal performance with large alignment datasets, mview uses threads and GPUs. 

#### Threads

*   **Linux**: OpenMP library (usually included with gcc as `libgomp`)
*   **macOS**: Install via Homebrew: `brew install libomp`
*   **Note**: alntools should compile without OpenMP

#### GPU Acceleration

On macOS with Apple Silicon, alntools uses GPUs to significantly boost performance.

*   **Automatic detection**: GPU acceleration is automatically enabled when Metal is available (run `system_profiler SPDisplaysDataType` to determine if you have GPUs)
*   **Performance**: 1000x faster bin queries for datasets with >10K alignments
*   **Fallback**: Automatically falls back to CPU if GPU not available

#### Install alntools

```bash
# Clone alntools into the same tools directory
cd /path/to/your/makeshift/tools
git clone https://github.com/eitanyaffe/alntools.git
cd alntools
make install
```

Add the local bin directory to your path if needed:
```
export PATH="$HOME/.local/bin:$PATH"
```

**Set MAKESHIFT_ROOT environment variable** (only if using alntools): mview needs to locate the alntools C++ source files to compile the R adapter. Set this to your makeshift directory:

```bash
export MAKESHIFT_ROOT=/path/to/your/makeshift
```

Add this to your shell profile (`.bashrc`, `.zshrc`, etc.) to make it permanent.

Install R package:
```r
install.packages("Rcpp")
```

**Note**: The alntools adapter for mview is compiled during the first run of mview.

## Quick Start

### 1. Start R and launch mview

```bash
# Navigate to the mview directory
cd /path/to/your/makeshift/tools/mview
# Start R
R
```

```r
# Launch mview
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
- **Bottom tabs**: Tables (regions, contigs, genomes, genes, alignments).

#### Selecting the viewing region

1. **Select contigs**: Use Regions, Contigs or Genomes tab to select one or more contigs
2. **Zoom**: Use mouse to select regions in the visualization
3. **Lock zoom**: Press Alt+Z to set the selected region as zoom
4. **Window size**: Press Alt+N (N=2-7) to zoom to 10^N basepairs around current center
5. **Navigation**: Use Alt+Backspace for undo, press Help button for all shortcuts

#### Tables

The bottom panel contains these tabs:
- **Regions**: Manage saved visualization regions 
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

### Profiles

mview includes several built-in profiles for genomic visualization:

- **alignment**: Display read alignments with advanced mutation analysis, multiple visualization modes, and comprehensive filtering options including alignment length, mutation rates, and clipping patterns
- **gene**: Display gene segments with strand orientation and customizable coloring
- **segments**: Display genomic segments as colored rectangles (used for regions, annotations)
- **axis**: Show contig coordinates with tick marks for orientation

Profiles support interactive parameter adjustment through the right panel and hover information for all elements. For detailed information about profile parameters and usage, see the comprehensive guide: [profiles](docs/profiles.md).

### Managing regions

The Regions tab provides comprehensive management of saved visualization states, including assembly, contig selection, and zoom level. Region tables are stored as tab-delimited files with automatic versioning.

#### Creating and opening region tables

- **New Table**: Create a new region table with a custom name
- **Open Table**: Select from existing region tables in the regions directory  
- **Delete Table**: Remove an entire region table (with confirmation)

#### Working with regions

**To save current region**:
1. Navigate to Regions tab
2. Click "Add" button or press Alt+S
3. Enter description in the dialog
4. Region is automatically saved

**To restore a region**:
1. Select region from the table
2. Click "Goto" button in the Actions column
3. Use Alt+Backspace to undo navigation

**To manage regions**:
- **Edit**: Modify region descriptions
- **Delete**: Remove individual regions

Region files are saved as `[name].txt` with automatic backup versions stored in `versions/[name]_v1.txt`, `v2.txt`, etc.

#### Regions profile visualization

Saved regions are automatically displayed in the segments profile, which shows each single-contig region as a colored rectangle with its ID label. This provides immediate visual feedback when navigating between regions and helps identify overlapping or adjacent saved areas. The regions profile appears between the gene and axis profiles in both the c60 and minimal configurations.

## Working with your data

See [data_setup.md](docs/data_setup.md) for detailed instructions on preparing data files and configuration.

## Customizing views and profiles

**Views**: Create custom views by copying existing view files and modifying which profiles to display and their parameters. Each view defines a specific layout of gene tracks, read alignments and other profiles.

**Profiles**: Building custom profiles is straightforward - use existing profile code from the `profiles/` directory as a template and modify the plotting functions to handle your specific data types and visualization needs. The segments profile has minimal code and can be used as a starting point, see [profiles](docs/profiles.md).

## Developement log

### Version v1.01

- Alignment bin mode enhanced with density of segregating sites, non-reference sites, and distribution of reads by logarithmic mutation distance categories.
- Enhanced alignment filtering with granular clipping modes and min/max mutations percentage range controls across all query modes.
- Alignment length filtering: Added minimum and maximum alignment length filters based on read coordinates for precise control over which alignments to include in visualization.
- GPU support added (Apple Silicon only) for bin queries in the alignment profile.
- Added nucleotide sequence display for axis profiles.

### Under development (v1.02)

- Legends tab added
- Short indel filtering: Added `min_indel_length` parameter to alignment profiles to filter out short indels from mutation density calculations, reducing noise from sequencing artifacts while preserving longer, more meaningful indels.
