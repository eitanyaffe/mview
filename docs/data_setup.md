# Working with your data

**Quick overview**: Create your data tables, copy the minimal configuration, edit paths and field selections.

## Step 1: Prepare your data tables

You need to create several tab-delimited text files that describe your assemblies and genomic features:

### Required files

**1. Assembly Table** (`assembly_table.txt`)
```
assembly_id
your_assembly_1
your_assembly_2
```

**2. Contig Tables** (`contig_table_{assembly_id}.txt`)
```
contig	length	circular
ctg001	50000	FALSE  
ctg002	75000	FALSE
ctg003	45000	TRUE
```
Fields: `contig` (ID), `length` (bp), `circular` (TRUE/FALSE)

### Required files (continued)

**3. Genome Tables** (`genome_table_{assembly_id}.txt`)
```
gid	length
genome_1	1500000
genome_2	2100000
```
Fields: `gid` (genome ID), `length` (total genome length)

**4. Contig Mapping Tables** (`contig_map_table_{assembly_id}.txt`)
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
rl("my_project")  # loads /path/to/your/configs/my_project/minimal_cfg.r
```

## Adding read alignment visualization

**Generate alignments**: 
```bash
minimap2 -x sr assembly.fasta reads_R1.fastq reads_R2.fastq > alignment.paf
```

**Convert to ALN format**:
```bash
alntools construct -contigs assembly.fasta -paf alignment.paf -output alignment.aln
```

**Configure view**: Add alignment profile to your view file pointing to ALN files. The ALN format supports multiple visualization modes based on zoom level.