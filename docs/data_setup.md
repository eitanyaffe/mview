# Working with your data

**Quick overview**: Create your data tables, copy the minimal configuration, edit paths and field selections.

## Step 1: Prepare your data tables

You need to create several tab-delimited text files that describe your assemblies and genomic features:

### Required files

**1. Assembly Table**
```
assembly_id
your_assembly_1
your_assembly_2
```

Each assembly requires an assembly table, a contig table and a contig mapping table.

**2. Contig Tables**


```
contig	length	circular
ctg001	50000	FALSE  
ctg002	75000	FALSE
ctg003	45000	TRUE
```
Fields: `contig` (ID), `length` (bp), `circular` (TRUE/FALSE)

**3. Genome Tables**
```
gid	length
genome_1	1500000
genome_2	2100000
```
Fields: `gid` (genome ID), `length` (total genome length)

**4. Contig Mapping Tables**
```
contig	gid
ctg001	genome_1
ctg002	genome_1  
ctg003	genome_2
```
Fields: `contig`, `gid` (links contigs to genomes)

**5. Gene Tables** (`gene_table.txt`, optional)
```
gene	assembly	contig	start	end	strand	prot_desc	tax	identity	coverage
gene001	BAA	ctg001	100	1200	+	DNA polymerase	Bacteria	85.5	95.2
```
**Mandatory fields**: `gene`, `assembly`, `contig`, `start`, `end`  
**Optional fields**: `strand`, `uniref`, `identity`, `coverage`, `evalue`, `bitscore`, `prot_desc`, `tax`, `uniref_count`

Edit the `keep` vector in your config to match your available fields.

### Creating alignment files

To visualize read alignments in mview, you need to align your reads to the assembly and convert the results to ALN format for efficient querying.

**Step 1: Create minimap2 index**
```bash
# Create index for faster alignment (recommended for large assemblies)
minimap2 -I 1G -x map-hifi \
        -t 8 \
        -d ref.mmi \
        contigs.fasta
```

**Step 2: Generate alignments**
```bash
# Align reads to assembly and generate PAF format
minimap2 -I 1G --cs -c -X \
        -N 10 \
        -t 8 \
        ref.mmi \
        reads.fastq \
        > alignment.paf
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
alntools construct -contigs contigs.fasta -paf alignment.paf -output alignment.aln
```

**Step 4: Configure view**: Add alignment profile to your view file pointing to ALN files. The ALN format supports multiple visualization modes (full, bin, pileup) based on zoom level.

## Step 2: Create your configuration

1. **Copy the minimal configuration**:
```bash
cp -r mview/configs/minimal /path/to/your/configs/my_project
```

2. **Edit the configuration file** `/path/to/your/configs/my_project/minimal_cfg.r`:
```r
# Update data directories
example_dir <- "/path/to/your/data"
tables_dir <- file.path(example_dir, "tables")

# Edit gene fields if needed
keep <- c("gene","contig","start","end","strand","prot_desc","tax")  # customize as needed
```

3. **Load your configuration**:
```r
source("mview.r")
rl("my_project", cdir="/path/to/your/configs")
```
