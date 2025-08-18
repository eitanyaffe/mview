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
install.packages(c("shiny", "DT", "plotly", "ggplot2", "shinyjs", "shinyjqui"))
```

### Installing alntools (optional; needed for read alignment profiles)

For optimal performance with large alignment datasets, mview uses threaded bin queries from alntools. This requires:

*   **Linux**: OpenMP library (usually included with gcc as `libgomp`)
*   **macOS**: Install via Homebrew: `brew install libomp`

#### GPU Acceleration (Apple Silicon)

On macOS with Apple Silicon, mview supports GPU acceleration for significant performance improvements with large alignment datasets:

*   **Automatic detection**: GPU acceleration is automatically enabled when Metal is available (run `system_profiler SPDisplaysDataType`)
*   **Performance**: 1000x faster bin queries for datasets with >10K alignments
*   **Fallback**: Automatically falls back to CPU if GPU not available

Note: alntools should compile and work without OpenMP, but large alignment datasets will process sequentially rather than in parallel.

If you plan to visualize read alignments, install alntools:

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
- **Bottom tabs**: Tables (states, contigs, genomes, genes, alignments).

#### Selecting the viewing region

1. **Select contigs**: Use State, Contigs or Genomes tab to select one or more contigs
2. **Zoom**: Use mouse to select regions in the visualization
3. **Lock zoom**: Press Alt+Z to set the selected region as zoom
4. **Window size**: Press Alt+N (N=2-7) to zoom to 10^N basepairs around current center
5. **Navigation**: Use Alt+Backspace for undo, press Help button for all shortcuts

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

### Profiles

- **alignment**: Display read alignments with multiple visualization modes:
  - **by_mut_density**: Traditional mutation density visualization (gray to red scale)
  - **by_seg_density**: Segregating sites density visualization (blue scale) 
  - **by_nonref_density**: Non-reference sites density visualization (yellow to orange scale)
  - **by_genomic_distance**: Stacked mutation distance categories (red gradient, log scale)
- **gene**: Display gene segments.

#### Alignment Profile Parameters

The alignment profile supports several key parameters accessible through the right panel:

- **bin_style**: Choose visualization mode (`by_mut_density`, `by_seg_density`, `by_nonref_density`, `by_genomic_distance`)
- **seg_threshold**: Threshold for segregating sites detection (default: 0.2)
- **non_ref_threshold**: Threshold for non-reference sites detection (default: 0.9)
- **plot_style**: Query mode selection (`auto_full`, `auto_pileup`, `bin`, `full`, `pileup`)
- **binsize/target_bins**: Control bin resolution for bin mode
- **use_gpu**: Enable GPU acceleration for bin queries (Apple Silicon only, default: true)

**Alignment Filtering**:
- **clip_mode**: Control which alignments to include based on clipping (`all`, `complete`, `allow_one_side_clip`, `only_one_side_clipped`, `only_two_side_clipped`)
- **clip_margin**: Margin of error when checking read boundaries (default: 10)
- **min_mutations_percent**: Minimum mutations percentage threshold with preset options (0%, 0.01%, 0.1%, 1%, 10%). Default 0%
- **max_mutations_percent**: Maximum mutations percentage threshold with preset options (0%, 0.01%, 0.1%, 1%, 10%). Default 10%

Hover shows element information in all profiles. Profile parameters can be set in code or adjusted interactively in the right panel. See the detailed profile guide: [profiles](docs/profiles.md).

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

## Current development version (v1.01)

- Alignment bin mode enhanced with density of segregating sites, non-reference sites, and distribution of reads by logarithmic mutation distance categories.
- Enhanced alignment filtering with granular clipping modes and min/max mutations percentage range controls across all query modes.