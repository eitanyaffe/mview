# Profiles

Profiles are the core visualization components in mview. Each profile displays a specific type of genomic data as a horizontal track that can be stacked with other profiles.

## Supported Profiles

### Gene Profile

Displays annotated genes as colored arrows or rectangles along the genome.

#### Key Parameters
- **`color_style`**: How to color genes
  - `"by_taxonomy"`: Color genes based on taxonomic classification
  - `"by_regex"`: Color genes based on regex pattern matching in descriptions
- **`show_labels`**: Whether to display gene names/IDs
- **`height_mode`**: How to size gene features
  - `"uniform"`: All genes same height
  - `"by_length"`: Gene height proportional to length
- **`strand_separation`**: Whether to separate forward/reverse strand genes vertically

#### Data Requirements
- Gene table with: `gene`, `contig`, `start`, `end`, `strand`
- Optional annotation table with: `gene`, `prot_desc`, `tax`

#### Example Usage
```r
gene_profile(
  id = "genes",
  name = "Genes", 
  height = 0.1,
  gene_f = get_genes_f,
  params = list(
    color_style = list(
      type = "select",
      choices = c("by_taxonomy", "by_regex"),
      default = "by_taxonomy"
    )
  )
)
```

### Alignment Profile

Shows read alignment data with coverage, variants, and alignment quality.

#### Key Parameters
- **`coverage_mode`**: How to display coverage
  - `"histogram"`: Coverage as histogram bars
  - `"line"`: Coverage as smooth line
  - `"both"`: Both histogram and line
- **`show_variants`**: Whether to highlight SNPs/indels
- **`quality_threshold`**: Minimum mapping quality to display
- **`max_reads`**: Maximum number of individual reads to show
- **`binning_size`**: Base pairs per coverage bin for performance

#### Data Requirements
- Alignment files (typically `.aln` format)
- Custom C++ alignment parsing (requires Rcpp)

#### Example Usage
```r
align_profile(
  id = "align_early",
  name = "Early Timepoint",
  aln_f = function(cxt) {
    tag <- get_map_tag(cxt$assembly, "early")
    get_data("MINIMAP_LIB_ALN", tag = tag, read_f = aln_load)
  },
  height = 0.4,
  params = list(
    coverage_mode = list(
      type = "select", 
      choices = c("histogram", "line", "both"),
      default = "histogram"
    ),
    show_variants = list(
      type = "boolean",
      default = TRUE
    )
  )
)
```

### Axis Profile

Provides coordinate reference scales and contig boundaries.

#### Key Parameters
- **`tick_interval`**: Spacing between coordinate ticks
- **`show_contig_names`**: Whether to label contig boundaries
- **`coordinate_format`**: How to format coordinate labels
  - `"bp"`: Raw base pairs (e.g., 12345)
  - `"kb"`: Kilobases (e.g., 12.3kb)
  - `"mb"`: Megabases (e.g., 1.2Mb)
- **`grid_lines`**: Whether to show vertical grid lines

#### Example Usage
```r
axis_profile(
  params = list(
    tick_interval = list(
      type = "integer",
      default = 10000
    ),
    coordinate_format = list(
      type = "select",
      choices = c("bp", "kb", "mb"),
      default = "kb"
    )
  )
)
```

## Profile Architecture

### Common Parameters
All profiles support these standard parameters:

- **`id`**: Unique identifier for the profile
- **`name`**: Display name shown in the interface
- **`height`**: Relative height (0.0-1.0) in the stacked view
- **`params`**: List of user-customizable parameters

### Parameter Types
- **`select`**: Dropdown menu with predefined choices
- **`boolean`**: On/off toggle
- **`integer`**: Numeric input for whole numbers
- **`double`**: Numeric input for decimal numbers
- **`string`**: Text input field

### Data Functions
Each profile typically requires a data function that:
1. Takes a context object (`cxt`) with assembly, contigs, and zoom info
2. Returns properly formatted data for visualization
3. Handles caching for performance

## Creating Custom Profiles

As an advanced feature, you can create your own profiles by:

1. **Following the profile API**: Implement required functions for data loading and plotting
2. **Defining parameters**: Specify user-controllable options
3. **Registering the profile**: Add to your view configuration

See [Advanced Topics](advanced.md#building-a-profile) for detailed instructions on creating custom profiles.

## Profile Performance

### Optimization Tips
- **Use caching**: Profile data functions should use the cache system
- **Limit data size**: Filter data to current zoom level when possible  
- **Bin large datasets**: Use binning for high-resolution data
- **Async loading**: Consider lazy loading for expensive computations

### Memory Management
- Profiles share the global cache system
- Data is automatically cached by genomic region
- Use `cache_info()` to monitor memory usage
- Clear cache between sessions if memory is limited