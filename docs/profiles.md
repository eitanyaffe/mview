## Profiles

This page describes the built-in profiles and their parameters.

### How to set profile parameters

- **In code (when defining views/profiles)**: pass constructor arguments and/or `params` to the profile factory functions (e.g., `align_profile(...)`, `gene_profile(...)`). Constructor arguments set defaults. The `params` list registers interactive controls in the right panel.
- **Interactively (in the app)**: use the right panel to change parameters at runtime. Values persist for the session (via cache), and updates are applied immediately to the plots.

Parameters listed in a profile's `params` appear in the right panel. Users can use this list to select which parameters are interactive.

Hover is available in all profiles to reveal information about the plotted elements.

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
  - `bin_type`: `auto` or a fixed bin size (e.g., `1000`).
  - `target_bins`: target number of bins when `bin_type=auto`.

#### Full mode parameters
  - `alignment_filter`:
    - `all`: no filtering.
    - `single`: include only reads with exactly one alignment in the queried intervals.
    - `single_complete`: like `single`, and the single alignment must span the entire read (no clipping).
    - `only_multiple`: include only reads with two or more alignments in the queried intervals.
  - `full_style`:
    - `none`: alignments gray; no mutation overlay.
    - `by_mutations`: color alignments by mutation density.
    - `by_strand`: color by strand orientation relative to the read.
    - `show_mutations`: alignments gray; draw per-mutation markers (limited by `max_mutations`).
  - `height_style`:
    - `by_coord_left`: minimize overlap by start coordinate.
    - `by_coord_right`: minimize overlap by end coordinate.
    - `by_mutations`: order by mutation density while preventing overlaps.
  - `max_reads`: cap alignments fetched per interval in full mode.
  - `max_mutations`: cap mutations drawn when `full_style=show_mutations`.
