## Profiles

This page describes the built-in profiles and their parameters.

### How to set profile parameters

- **In code (when defining views/profiles)**: pass constructor arguments and/or `params` to the profile factory functions (e.g., `align_profile(...)`, `gene_profile(...)`). Constructor arguments set defaults. The `params` list registers interactive controls in the right panel.
- **Interactively (in the app)**: use the right panel to change parameters at runtime. Values persist for the session (via cache), and updates are applied immediately to the plots.

Which parameters appear interactively is up to you: only parameters listed in a profile's `params` are exposed in the right panel. Include just the controls you want users to adjust for convenience; other parameters can remain constructor-only defaults.

Hover is available in profiles to reveal information about the plotted elements.

---

### alignment profile

The alignment profile adapts to zoom level and can also be forced to a specific mode.

What is shown:
- **bin mode**: per-bin coverage bars; fill color indicates mutation density (gray to red). Binning adapts to zoom unless a fixed bin size is selected.
- **full mode**: read spans and alignment rectangles; optional mutation overlay when enabled; alignment coloring indicates mutation density.
- **pileup mode**: stacked bars at each position showing per-variant counts; consistent variant colors across modes.

Parameter details
- **alignment_filter**:
  - `all`: no filtering.
  - `single`: include only reads that have exactly one alignment within the queried intervals.
  - `single_complete`: like `single`, but the single alignment must span the entire read (no clipping; read_start = 0 and read_end = read_length).
  - `only_multiple`: include only reads that have two or more alignments within the queried intervals.
- **height_style**:
  - `by_coord_left`: assign vertical lanes to minimize overlap by sorting reads by start coordinate.
  - `by_coord_right`: assign vertical lanes to minimize overlap by sorting reads by end coordinate.
  - `by_mutations`: sort reads by mutation density (highest first) and assign lanes while preventing overlaps.
- **max_reads**: applied at query time as a cap on the number of alignments fetched per interval for full mode, effectively limiting the number of reads rendered.
Main parameters
- **plot_style**: display strategy: `auto_full`, `auto_pileup`, `bin`, `full`, `pileup`.
- **alignment_filter**: subset alignments: `all`, `single`, `single_complete`, `only_multiple`.
- **full_style**: coloring/overlay in full mode: `none`, `by_mutations`, `by_strand`, `show_mutations`.
- **height_style**: how vertical stacking is determined in full mode: `by_mutations`, `by_coord_left`, `by_coord_right`.
- **target_bins**: target number of bins for automatic bin sizing in bin mode (integer).
- **full_threshold**: max visible range (bp) to use full mode with `auto_full` (integer bp).
- **pileup_threshold**: max visible range (bp) to use pileup mode with `auto_pileup` (integer bp).
- **max_reads**: cap on reads queried/rendered in full mode (integer).
- **max_mutations**: cap on mutations displayed when `show_mutations` is used (integer).
- **bin_type**: bin size selection in bin mode: `auto`, or fixed sizes like `10`, `100`, `1000`, `5000`, `10000`.

Notes
- Auto modes choose between bin/full or bin/pileup based on the current visible range and the corresponding thresholds.
- Variant colors are consistent across full and pileup modes.

---

### gene profile

What is shown:
- **simple mode** (zoomed out): compact markers at gene starts for efficiency.
- **full mode** (zoomed in): gene rectangles clipped to view; a short vertical segment at the promoter side indicates strand; coloring taken from a chosen column.

Main parameters
- **threshold**: range (bp) cutoff for switching simple vs full (integer bp).
- **color_field**: column name that determines the color (e.g., `tax_color`). If the provided name lacks `_color`, the suffix is appended.
- **height**: vertical size of the profile (pixels).

Notes
- If the color field is missing or all-NA, a default light gray is used.
- When many genes are visible in simple mode, sampling may be applied for responsiveness.

---

### axis profile

Shows the contig name, axis line, and tick marks for orientation.


