## Profiles

This page describes the built-in profiles and their parameters.

### How to set profile parameters

- **In code (when defining views/profiles)**: pass constructor arguments and/or `params` to the profile factory functions (e.g., `align_profile(...)`, `gene_profile(...)`). Constructor arguments set defaults. The `params` list registers interactive controls in the right panel.
- **Interactively (in the app)**: use the right panel to change parameters at runtime. Values persist for the session (via cache), and updates are applied immediately to the plots.

Parameters listed in a profile's `params` appear in the right panel. Users can use this list to select which parameters are interactive.

Hover is available in all profiles to reveal information about the plotted elements.

---

### segments profile

Display genomic segments as colored rectangles with customizable data sources. Designed for visualizing regions, annotations, or any segmented genomic features.

#### Main parameters
- **segments_f**: data source, can be:
  - Function: takes assembly as parameter and returns segment data (e.g., `function(assembly) get_data("CAV_REFINE_SEGMENT_TABLE", tag = assembly)`)
  - Cache key (string): references cached segment data (e.g., `"segments.current_regions"`)
  - Data frame: direct segment data with required columns (`assembly`, `contig`, `start`, `end`, `desc`, `id`)
- **color**: fill color for segment rectangles (default: `"#2E86AB"`)
- **height**: vertical size of the profile in pixels (default: 40)

#### Data format
Segment data requires these columns:
- **assembly**: assembly identifier (e.g., "EBC", "BAA")  
- **contig**: contig name within the assembly
- **start**: start coordinate within contig (local coordinates)
- **end**: end coordinate within contig (local coordinates)
- **desc**: description text for hover display
- **id**: unique identifier displayed as label

#### Display behavior
- **Assembly filtering**: only shows segments matching current assembly context
- **Context filtering**: automatically filters segments to current view and contigs
- **Rectangle display**: segments shown as filled rectangles with black borders
- **ID labels**: segment IDs displayed to the right of each rectangle in black text
- **Hover information**: shows contig coordinates and description when hovering over ID labels

#### Usage examples
```r
# Assembly segments with function
segments_profile(
  id = "assembly_segments",
  name = "Assembly Segments",
  segments_f = function(assembly) {
    get_data("CAV_REFINE_SEGMENT_TABLE", tag = assembly)
  }
)

# Basic segments profile with cache data
segments_profile(
  id = "regions",
  name = "Regions", 
  segments_f = "segments.current_regions"
)

# Custom segments with direct data and styling
segments_profile(
  id = "annotations",
  name = "Annotations",
  segments_f = my_segments_df,
  color = "#FF6B6B",
  height = 30
)
```

**Notes**:
- Supports string-based IDs for flexible labeling (e.g., "A1", "control", "REG_001")
- Automatically filters to show only segments relevant to the current assembly and view

---

### axis profile

Shows the contig name, axis line, and tick marks for orientation.

---

### gene profile

What is shown:
- **simple mode** (zoomed out): compact markers at gene starts for efficiency.
- **full mode** (zoomed in): gene rectangles clipped to view; a short vertical segment at the promoter side indicates strand; coloring taken from a chosen column.

Main parameters
- **threshold**: range (bp) cutoff for switching simple vs full (integer bp).
- **color_field**: column name that determines the color (e.g., `tax_color`). If the provided name lacks `_color`, the suffix is appended.
- **label_field**: column used for tooltip text. If missing or not provided, a default label is used.
- **height**: vertical size of the profile (pixels).

Notes
- When many genes are visible in simple mode, sampling may be applied for responsiveness.

**Coloring by multiple fields**. Create a color column from multiple gene fields (e.g., combine fields and map to hex colors), then set `color_field` to that column name.

---

### rrna profile

Display rRNA genes from barrnap GFF output with automatic type-based coloring and strand visualization.

#### Data requirements
The `get_gff_f` function should return a data frame with parsed GFF columns:
- **contig**: contig/sequence name
- **start/end**: gene coordinates  
- **strand**: strand direction (+ or -)
- **name**: rRNA gene name (e.g., "16S_rRNA")
- **product**: gene product description
- **rrna_type**: rRNA type for coloring (e.g., "16S", "23S", "5S")

### alignment profile

Plot read alignments. Depending on zoom or explicit settings, it renders coverage bins, full read alignments, or per-position pileups.

#### General parameters
- **plot_style**
  - `bin`: aggregate coverage into bins and draw coverage bars.
  - `full`: draw per-read alignments; can color by attributes or overlay per-read mutations.
  - `pileup`: draw stacked per-variant counts at each genomic position.
  - `auto_full`: use full when visible range ≤ `full_threshold`, otherwise bin.
  - `auto_pileup`: use pileup when visible range ≤ `pileup_threshold`, otherwise bin.
- **full_threshold**: bp cutoff for `auto_full`.
- **pileup_threshold**: bp cutoff for `auto_pileup`.

#### Bin mode parameters
  - `bin_type`: `auto` or a fixed bin size from predefined options (`1000`, `2000`, `5000`, `10000`, `20000`, `50000`, `100000`, `200000`, `500000`).
  - `target_bins`: target number of bins when `bin_type=auto`. Auto mode selects the largest available binsize ≤ `range_bp / target_bins`.
  - `bin_style`: visualization mode for bins:
    - `by_mut_density`: traditional mutation density (gray to red scale).
    - `by_median_mutation_density`: median mutation density across alignments per bin (gray to red scale).
    - `by_seg_density`: segregating sites density (blue scale).
    - `by_nonref_density`: non-reference sites density (yellow to orange scale).
    - `by_seg_clip_density`: segregating clip sites density (purple scale).
    - `by_non_ref_clip_density`: non-reference clip sites density (green scale).
    - `by_genomic_distance`: stacked mutation distance categories (red gradient, logarithmic scale).
  - `seg_threshold`: threshold for segregating sites detection (default: 0.2).
  - `non_ref_threshold`: threshold for non-reference sites detection (default: 0.9).


#### Mutation distance categories (logarithmic scale)
When using `bin_style=by_genomic_distance`, alignments are categorized by mutations per bp:
  - `dist_none`: exactly 0 mutations
  - `dist_5`: 1e-5 to 1e-4 per bp (10-100 mutations per 100kb)
  - `dist_4`: 1e-4 to 1e-3 per bp (100-1,000 mutations per 100kb)
  - `dist_3`: 1e-3 to 1e-2 per bp (1,000-10,000 mutations per 100kb)
  - `dist_2`: 1e-2 to 1e-1 per bp (10,000-100,000 mutations per 100kb)
  - `dist_1_plus`: >1e-1 per bp (>100,000 mutations per 100kb)

Categories are stacked from highest distance (bottom) to lowest distance (top) in shades of red.

#### Alignment filtering parameters
  - `clip_mode`: Controls which alignments are included based on clipping:
    - `all`: Allow all alignments regardless of clipping (default).
    - `complete`: Alignment must cover the entire read from start to end.
    - `allow_one_side_clip`: Allow alignments clipped on one side (start at read start OR end at read end).
    - `only_one_side_clipped`: Only include alignments clipped on exactly one side.
    - `only_two_side_clipped`: Only include alignments clipped on both sides.
    - `only_clipped`: Only include alignments clipped on one or both sides.
  - `clip_margin`: Margin of error in bases when checking read start/end positions (default 10).
  - `min_mutations_percent`: Minimum mutations percentage threshold with preset options (0%, 0.01%, 0.1%, 1%, 10%). Default 0%. Filters out alignments with mutation rate below this threshold.
  - `max_mutations_percent`: Maximum mutations percentage threshold with preset options (0%, 0.01%, 0.1%, 1%, 10%). Default 10%. Special case: 0% = only alignments with zero mutations. Otherwise filters out alignments with mutation rate above this threshold.
  - `min_indel_length`: Minimum indel length to include in mutation density calculations (default 3). Indels shorter than this threshold are filtered out from mutation counts and density calculations to reduce noise from sequencing artifacts. Set to 0 to include all indels.

#### Chunk-based height calculation

In full mode, alignments are grouped into "chunks" (stretches of related alignments) for height calculation, providing better visual separation while maintaining biological context. Heights are calculated for chunks, then inherited by individual alignments within each chunk. Clicking on any alignment still navigates to the full read context.

#### Full mode parameters
  - `full_style`:
    - `none`: alignments gray; no mutation overlay.
    - `by_mutations`: color alignments by mutation density.
    - `by_strand`: color by strand orientation relative to the read.
    - `show_mutations`: alignments gray; draw per-mutation markers (limited by `max_mutations`).
  - `height_style`:
    - `by_coord_left`: minimize overlap by start coordinate.
    - `by_coord_right`: minimize overlap by end coordinate.
    - `by_mutations`: order by mutation density while preventing overlaps.
  - `chunk_type`: defines how alignments are grouped into chunks for height calculation:
    - `break_on_overlap`: start new chunk when alignments overlap in read coordinates (default).
    - `break_on_gap`: start new chunk when gap between alignments in read coordinates exceeds max_gap.
    - `read`: entire read forms one chunk (equivalent to read-based heights).
    - `alignment`: each alignment forms its own chunk (maximum granularity).
  - `max_gap`: maximum gap tolerance for chunk detection in read coordinates (default: 10).
  - `max_reads`: cap alignments fetched per interval in full mode.
  - `max_mutations`: cap mutations drawn when `full_style=show_mutations`.

#### Mutation Identifiers

Mutation identifiers follow standard notation:
- **X:Y** - substitution where base X is substituted by base Y
- **+XYZ** - insertion where sequence XYZ was added
- **-XYZ** - deletion where sequence XYZ was removed

---

### synteny profile

Visualizes sequencing coverage and mutation density across genomic bins for multiple libraries. Shows which libraries have sufficient coverage in each genomic region and can display mutation patterns.

#### Main parameters
- **style**: visualization mode
  - `summary`: shows count of qualifying libraries per genomic bin (bar chart)
  - `detail`: shows heatmap with libraries as rows and genomic bins as columns
- **binsize**: bin size for analysis
  - `auto`: automatically selects optimal binsize based on view range and target_bins
  - Fixed sizes: `1000`, `2000`, `5000`, `10000`, `20000`, `50000`, `100000`, `200000`, `500000`
- **target_bins**: target number of bins when binsize is `auto` (default: 200)
- **min_xcov**: minimum x-coverage threshold (sequenced_bp / binsize) for inclusion (default: 1.0)
- **color_style**: coloring scheme
  - `none`: gray coloring for all qualifying data points
  - `mutations`: colored by mutation density using same red scale as alignment profile
- **hide_self**: whether to hide libraries that start with current assembly ID (default: true)
- **height**: vertical size of the profile in pixels (default: 400)

#### Data requirements
The synteny profile requires two input matrices with identical structure:
1. **Sequenced bp matrix**: sequencing coverage data (bp sequenced per bin per library)
2. **Mutation density matrix**: median mutation density per bin per library

Both matrices must have:
- First 3 columns: `contig`, `start`, `end` (genomic coordinates)
- Remaining columns: library data (one column per library)

#### Data function
The profile expects a simple `synteny_f` function that loads individual tables:
```r
get_synteny_f <- function(assembly, field, binsize, hide_self = TRUE) {
  data <- get_data("MINIMAP_SYNTENY_ASSEMBLY_TABLE", 
                   tag = paste0(assembly, "_", field, "_", binsize),
                   null.on.missing = TRUE)
  
  # filter out self-assembly libraries if hide_self is TRUE and data exists
  if (hide_self && !is.null(data)) {
    keep_cols <- !grepl(paste0("^", assembly, "_"), colnames(data))
    keep_cols[1:3] <- TRUE  # always keep contig, start, end
    data <- data[, keep_cols, drop = FALSE]
  }
  
  return(data)
}
```

The profile automatically:
- Calls this function twice: once for `"sequenced_bp"` and once for `"median_mutation_density"`
- Passes the `hide_self` parameter to the user function
- Handles missing data gracefully (returns `NULL` if either dataset is missing)
- Validates data structure and provides informative error messages

#### Display behavior
- **Summary mode**: Y-axis shows count of libraries meeting min_xcov threshold per bin
- **Detail mode**: Y-axis shows library names, X-axis shows genomic position; only displays bins where library meets threshold
- **X-coverage calculation**: `sequenced_bp / binsize` determines if library qualifies for each bin
- **Automatic binsize selection**: chooses largest available binsize ≤ `range_bp / target_bins`
- **Hover information**: shows position, coverage, mutation density, and library details

#### Usage example
```r
# Load synteny profile components
source("profiles/synteny/synteny_profile.r")
source("profiles/synteny/synteny_profile_detail.r") 
source("profiles/synteny/synteny_profile_summary.r")

# Create synteny profile
synteny_profile(
  id = "synteny",
  name = "synteny",
  synteny_f = get_synteny_f,
  params = default_synteny_params
)
```

**Notes**:
- Uses standard mview genomic coordinate mapping and filtering
- Integrates with parameter system for interactive control
- Compatible with mview's caching and view management