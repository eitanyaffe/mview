# Working with your data

**Quick overview**: Create your data tables, copy the minimal configuration, edit paths and field selections.

> **Note:** The column names and file structures shown below are examples only. Your existing data can have different column names and structures. You can use your existing tables through configurable "get functions" that transform tables as needed. Read through the entire guide before creating new files to understand how the system works.

## Segments: the core unit of visualization

mview uses a **segment-based** coordinate system. A segment is a sub-contig interval defined by a contig name, start, and end coordinates. A genome (or bin) is a collection of segments. The current view is an ordered sequence of segments that determines what is displayed.

You have two options for defining segments:

**Option 1: Simple 1-to-1 mapping** (recommended for most users)  
Define only a contig table. Your configuration will auto-generate segments where each segment spans an entire contig. This is the approach used in the minimal example.

**Option 2: Custom segment table**  
Provide your own segment table with sub-contig intervals. This is useful when you want to:
- Display arbitrary sub-regions of contigs
- Segment contigs based on alignment breakpoints
- Reorder segments independently of genomic coordinates

If you want to segment your contigs, alntools provides segmentation based on partially aligned reads (see [alntools documentation](https://github.com/eitanyaffe/alntools)).

## Step 1: Prepare your data tables

You need to have several tab-delimited text files that describe your assemblies and genomic features. The examples below show the required information, but your actual files can have different column names and structures:

**1. Assembly Table**
```
assembly_id
your_assembly_1
your_assembly_2
```

Each assembly requires an assembly table, a contig table, and a genome table. For simple 1-to-1 segments, add a contig mapping table. For custom segments, add a segment table and segment mapping table.

**2. Contig Tables**

```
contig	length	circular
ctg001	50000	FALSE  
ctg002	75000	FALSE
ctg003	45000	TRUE
```
Required information: contig identifier, sequence length. The circular field is optional.

**3. Genome Tables**
```
gid	length
genome_1	1500000
genome_2	2100000
```
Required information: genome identifier, total genome length

**4. Contig Mapping Tables** (for simple 1-to-1 segments)
```
contig	gid
ctg001	genome_1
ctg002	genome_1  
ctg003	genome_2
```
Required information: contig identifier, genome identifier (links contigs to genomes)

When using simple 1-to-1 mapping, your configuration will convert this to a segment map automatically.

**5. Segment Tables** (optional, for custom segments)
```
segment	contig	start	end	length
seg001	ctg001	1	25000	25000
seg002	ctg001	25001	50000	24999
seg003	ctg002	1	75000	75000
```
Required information: segment identifier, contig, start position, end position, length

**6. Segment Mapping Tables** (required if using custom segments)
```
segment	gid
seg001	genome_1
seg002	genome_1
seg003	genome_2
```
Required information: segment identifier, genome identifier (links segments to genomes)

**7. Gene Tables**
```
gene	assembly	contig	start	end	strand	prot_desc	tax	identity	coverage
gene001	BAA	ctg001	100	1200	+	DNA polymerase	Bacteria	85.5	95.2
```
**Required information**: gene identifier, assembly, contig, start position, end position  
**Optional information**: strand, protein description, taxonomy, identity, coverage, etc.

### Creating alignment files

To visualize read alignments in mview, you need to align your reads to the assembly and convert the results to ALN format for efficient querying. Alntools is a C++ program that converts alignment files (PAF format) into a binary format optimized for real-time visualization ([github.com/eitanyaffe/alntools](https://github.com/eitanyaffe/alntools)). The R interface provides direct access to memory-stored alignments, enabling rapid interactive exploration of millions of reads without performance bottlenecks.

**Step 1: Create minimap2 index**
```bash
# Create index for faster alignment (recommended for large assemblies)
minimap2 -I 1G -x map-hifi -t 8 -d ref.mmi contigs.fasta
```

**Step 2: Generate alignments**
```bash
# Align reads to assembly and generate PAF format
minimap2 -I 1G --cs -c -X -N 10 -t 8 ref.mmi reads.fastq > alignment.paf
```

Parameters explained:
- `-I 1G`: Load at most 1GB of reference sequences into memory
- `--cs`: Generate CIGAR strings for detailed alignment information
- `-c`: Generate CIGAR in PAF output
- `-X`: Include sequence in PAF for better visualization
- `-N 10`: Find up to 10 secondary alignments
- `-t 8`: Use 8 threads

**Step 3: Convert to ALN format**
```bash
# Convert PAF to ALN format for efficient querying in mview
alntools construct -ifn_paf alignment.paf -ofn alignment.aln
```

## Step 2: Create your configuration

### Setting up your configuration

1. **Create project directory and copy configuration files**:
```bash
# Create your project directory
mkdir -p /path/to/your/configs/my_project

# Copy minimal configuration files
cp configs/minimal/minimal_cfg.r /path/to/your/configs/my_project/my_project_cfg.r
cp configs/minimal/minimal_view.r /path/to/your/configs/my_project/my_project_view.r
```

2. **Edit the configuration file** `/path/to/your/configs/my_project/my_project_cfg.r`:
```r
# Update data directories to point to your data
example_dir <- "/path/to/your/data"
tables_dir <- file.path(example_dir, "tables")

# Modify get functions to match your data structure
get_assembly_f <- function() {
  # Point to your assembly table
  read_cached("assembly_table", "/your/actual/path/assemblies.txt")
}

# Edit other get functions similarly, see more on get functions below

# Update view registration to point to your view file
view_file <- "/path/to/your/configs/my_project/my_project_view.r"
view_register("example", view_file)
```

3. **Load your configuration**:
```r
source("mview.r")
rl("my_project", cdir="/path/to/your/configs")
```

### Understanding get functions

The configuration file uses "get functions" that tell mview how to access your data. Your data can be stored anywhere and in any format. The get functions act as adapters between your actual data files and mview's expected format. You can transform field names, filter data, or compute required fields on the fly.

**Example**: If your assembly table has a field called `subjects` instead of `assembly_id`, your get function can transform it:

```r
get_assembly_f <- function() {
  df <- read.delim("/path/to/your/data.txt")
  # Transform field name
  data.frame(assembly_id = df$subjects)
}
```

### Required get functions

Each configuration must define these core get functions:

**1. `get_assembly_f()`** - Returns available assemblies
```r
get_assembly_f <- function() {
  # Must return: assembly_id column
  read.delim("/your/path/assembly_table.txt")
}
```

**2. `get_contigs_f(assembly)`** - Returns contigs for an assembly
```r
get_contigs_f <- function(assembly = NULL) {
  # Must return: contig, length, circular columns
  path <- paste0("/your/path/contigs_", assembly, ".txt")
  df <- read.delim(path)
  # Transform if needed - example renaming 'size' to 'length'
  data.frame(contig = df$contig_id, length = df$size, circular = df$is_circular)
}
```

**3. `get_genomes_f(assembly)`** - Returns reconstructed genomes
```r
get_genomes_f <- function(assembly = NULL) {
  # Must return: gid, length columns
  df <- read.delim(paste0("/your/path/genomes_", assembly, ".txt"))
  # Can compute length on the fly
  df$length <- df$length_kb * 1000  # example conversion
  df[, c("gid", "length")]
}
```

**4. `get_segments_f(assembly)`** - Returns segments for an assembly

For simple 1-to-1 mapping (auto-generate from contigs):
```r
get_segments_f <- function(assembly = NULL) {
  # Must return: segment, contig, start, end, length columns
  contigs <- get_contigs_f(assembly)
  if (is.null(contigs)) return(NULL)
  data.frame(
    segment = paste0("s", seq_len(nrow(contigs))),
    contig = contigs$contig,
    start = 1L,
    end = contigs$length,
    length = contigs$length,
    stringsAsFactors = FALSE
  )
}
```

For custom segment tables:
```r
get_segments_f <- function(assembly = NULL) {
  # Must return: segment, contig, start, end, length columns
  read.delim(paste0("/your/path/segments_", assembly, ".txt"))
}
```

**5. `get_segment_map_f(assembly)`** - Maps segments to genomes

For simple 1-to-1 mapping (auto-generate from contig map):
```r
get_segment_map_f <- function(assembly = NULL) {
  # Must return: segment, gid columns
  segments <- get_segments_f(assembly)
  if (is.null(segments)) return(NULL)
  
  # read contig->gid mapping
  contig_map <- read.delim(paste0("/your/path/contig_map_", assembly, ".txt"))
  
  # convert to segment->gid mapping
  ix <- match(contig_map$contig, segments$contig)
  data.frame(segment = segments$segment[ix], gid = contig_map$gid, stringsAsFactors = FALSE)
}
```

For custom segment tables:
```r
get_segment_map_f <- function(assembly = NULL) {
  # Must return: segment, gid columns
  read.delim(paste0("/your/path/segment_map_", assembly, ".txt"))
}
```

### Optional get functions

**6. `get_genes_f(cxt)`** - Returns gene annotations (for gene tab)
```r
get_genes_f <- function(cxt) {
  # Must return: gene, contig, start, end columns
  # Optional: strand, prot_desc, tax, identity, coverage, etc.
  genes <- read.delim("/your/path/genes.txt")
  genes <- genes[genes$assembly_name == cxt$assembly, ]  # filter by assembly
  
  # Transform field names as needed
  genes$gene <- genes$gene_id
  genes$prot_desc <- genes$annotation
  
  # Select only needed fields
  fields <- c("gene", "contig", "start", "end", "strand", "prot_desc", "tax")
  genes[, intersect(names(genes), fields)]
}
```

**7. `get_aln_f(assembly)`** - Returns alignment data (for alignment profiles)
```r
get_aln_f <- function(assembly = NULL) {
  # Returns alignment object for visualization
  path <- paste0("/your/path/alignments_", assembly, ".aln")
  aln_load(path)  # uses alntools to load ALN format
}
```
