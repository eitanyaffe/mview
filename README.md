# Mview

A nimble interactive genomic data visualization tool built with R Shiny for exploring assemblies, contigs, genes, and genomic features.

## Description

Mview provides a highly customizable R-based platform for visualizing genomic data through user-defined interactive views. Mview allows to explore read alignments and gene annotations. The tool is designed for researchers working with reference-based and reconstructed metagenomic assemblies.

## Prerequisites

### Installation

Install all required packages at once:

```r
install.packages(c("shiny", "DT", "plotly", "ggplot2", "shinyjs", "shinyjqui"))

# Required to visualize read alignments using alntools
install.packages(c("Rcpp"))
```

## Quick Start

### 1. Start R and Launch Mview

```r
# Start R in the mview directory
source("app.r")
```

### 2. Load a Configuration

Once the app starts, you need to load a configuration. In the R console, type:

```r
# Load the simple example configuration
rl("example_simple")

# Or load your own custom configuration
rl("my_config")
```

### 3. Use the Application

The app will open in your browser. You can then:
- Select assemblies and contigs from the data tables
- Explore gene annotations through interactive profiles
- Use interactive plots to zoom and navigate genomic regions
- Use keyboard shortcuts ('r' to reset zoom, 'q' to quit)
- Adjust visualization parameters in the collapsible parameter panel

### What `rl()` Does

The `rl()` function (reload) loads a configuration by:
1. Clearing any existing cache and state
2. Sourcing the configuration file from `configs/{config_id}/{config_id}_cfg.r`
3. Registering all profiles, views, and tabs
4. Updating the UI with the new configuration

## UI Parameters

mview provides interactive parameters to customize how profiles are displayed. These parameters appear in a collapsible panel on the right side of the interface.

### Using Parameters
- **Parameter Panel**: Click the toggle button to expand/collapse the parameter panel
- **Profile-Specific Settings**: Each profile can define its own customizable parameters
- **Real-Time Updates**: Changes to parameters immediately update the visualization
- **Parameter Types**: 
  - Select dropdowns (e.g., color schemes)
  - Numeric inputs (e.g., thresholds, sizes)
  - Boolean toggles (e.g., show/hide features)

### Example Parameters
- **Gene Profile**: Color genes by taxonomy vs. regex patterns
- **Custom Profiles**: Any parameter you define in your profile configuration

## Saving and Loading States

mview automatically manages session state and provides persistence across sessions.

### Automatic State Management
- **Current Selection**: Assembly and contig selections are maintained
- **Zoom Level**: Current genomic coordinate range is preserved
- **Parameter Values**: All parameter settings are cached
- **View State**: Active view and profile configurations are remembered

### State Persistence
- **Session Cache**: State is maintained within an R session
- **Cross-Session**: Parameter values persist between app restarts through the cache system
- **Configuration-Specific**: Each configuration maintains its own cached state

### Manual State Control
- **Reset**: Use 'r' keyboard shortcut to reset zoom to show all contigs
- **Clear Cache**: Restart R session to completely clear all cached state
- **Configuration Switch**: Use `rl("config_id")` to switch between different configurations

## How Mview is Structured

### State Management
mview maintains application state through reactive variables:
- **Assembly**: Currently selected assembly ID
- **Contigs**: Selected contigs for visualization
- **Zoom**: Current genomic coordinate range for detailed views
- **Parameters**: Profile-specific display options

### Profiles
Profiles are the core visualization components that stack vertically:
- **Gene Profile**: Displays annotated genes with taxonomy coloring
- **Alignment Profile**: Shows read alignments with coverage and variants
- **Axis Profile**: Provides coordinate reference scales
- **Custom Profiles**: Advanced users can create profiles for any genomic data type

Each profile includes customizable parameters for colors, filtering, and display options.

📖 **See [Available Profiles](docs/profiles.md) for detailed documentation and parameters**

### Tabs
Tabs provide data exploration interfaces:
- **Genes Tab**: Browse and filter gene annotations
- **Alignments Tab**: Explore read mapping data (when available)

Data tables are interactive with sorting, filtering, and selection capabilities.

## Keyboard Shortcuts

- **`r`**: Reset zoom to show all selected contigs
- **`q`**: Quit application
- **Arrow keys**: Navigate within tables and interface elements
- **Mouse wheel**: Zoom in/out on plots
- **Click + drag**: Pan across genomic regions

## Configuration Files

mview uses a modular configuration system:

```
configs/
├── example_simple/           # Simple example for learning
│   ├── example_simple_cfg.r  # Main configuration
│   └── example_simple_main_view.r  # Profile definitions
└── my_config/                # Your custom configuration
    ├── my_config_cfg.r       # Your configuration
    └── my_config_main_view.r # Your profiles
```

## Data Structure

```
examples/
├── lookup.txt                # Data table registry
└── tables/                   # Data files
    ├── assembly_table.txt    # Assembly definitions
    ├── contig_table_*.txt    # Contig information per assembly
    ├── gene_table_*.txt      # Gene annotations per assembly
    └── uniref_table_*.txt    # Protein annotations per assembly
```

## Advanced Topics

For detailed information on customization and advanced usage, see:

### [Understanding Views and Profiles](docs/advanced.md#understanding-views-and-profiles)
Learn how mview's visualization system works with views and profiles.

### [Importing Your Own Data](docs/advanced.md#importing-your-own-data)  
Step-by-step guide to import and visualize your genomic data.

### [Available Profiles](docs/profiles.md)
Complete documentation of current profiles and their parameters.

### [Building Custom Profiles](docs/advanced.md#building-a-profile)
Advanced: Create custom visualization profiles for your specific data types.

### [Data Requirements](docs/advanced.md#data-requirements)
Required table formats and column specifications.

### [Configuration Architecture](docs/advanced.md#configuration-architecture)
Deep dive into mview's modular architecture and customization points.