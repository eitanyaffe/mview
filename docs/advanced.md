# Advanced Topics

## Understanding Views and Profiles

Mview's visualization system is built around two key concepts:

### Views
- **Definition**: A view defines a complete visualization setup for a specific analysis context
- **Contents**: Collections of related profiles that work together
- **Selection**: Users switch between views via the interface

### Profiles  
- **Definition**: Individual visualization tracks that display specific types of genomic data
- **Stacking**: Multiple profiles stack vertically to create comprehensive views
- **Examples**:
  - `gene_profile()`: Displays annotated genes with customizable coloring
  - `align_profile()`: Shows read alignment data with coverage and variants
  - `axis_profile()`: Provides coordinate reference scales
- **Customization**: Each profile can define its own parameters for user control

### Relationship
```
Configuration
├── View 1 (e.g., "Early timepoint")
│   ├── Gene Profile
│   ├── Alignment Profile  
│   └── Axis Profile
└── View 2 (e.g., "Late timepoint")
    ├── Gene Profile
    ├── Different Alignment Profile
    └── Axis Profile
```

## Importing Your Own Data

To use mview with your own genomic data, follow these steps:

### 1. Start with the Example
```bash
cp -r configs/example_simple configs/my_config
```

### 2. Prepare Your Data Tables
Create your data tables anywhere on your system. Required tables:
- **Assembly table**: Lists your assemblies/samples
- **Contig tables**: Contig information for each assembly  
- **Gene tables**: Gene annotations for each assembly
- **Annotation tables**: Protein/functional annotations (optional)

### 3. Create a Lookup Table
Create a single lookup file pointing to all your data:
```
# my_data/lookup.txt
id	path
ASSEMBLY_TABLE	/path/to/my_assemblies.txt
ASSEMBLY_CONTIG_TABLE:sample1	/path/to/sample1_contigs.txt
ASSEMBLY_CONTIG_TABLE:sample2	/path/to/sample2_contigs.txt
PRODIGAL_GENE_TABLE:sample1	/path/to/sample1_genes.txt
PRODIGAL_GENE_TABLE:sample2	/path/to/sample2_genes.txt
```

### 4. Customize Configuration Files
Edit `configs/my_config/my_config_cfg.r`:
- Update `set_lookup("path/to/your/lookup.txt")`
- Modify data loading functions to match your table column names
- Adjust assembly and contig registration functions

Edit `configs/my_config/my_config_main_view.r`:
- Customize profiles for your specific data types
- Adjust gene loading functions to match your annotation fields
- Add or remove profiles as needed

### 5. Load Your Configuration
```r
rl("my_config")
```

See `configs/example_simple/` for a complete working example with all required files and functions.

## Building a Profile

Profiles define how genomic features are visualized:

```r
# Example gene profile
gene_profile(
  id = "my_genes",              # Unique identifier
  name = "My Genes",            # Display name
  height = 0.1,                 # Relative height (0.0-1.0)
  gene_f = my_gene_function,    # Function returning gene data
  params = list(               # Customizable parameters
    color_style = list(
      type = "select",
      choices = c("by_taxonomy", "by_function"),
      default = "by_taxonomy"
    )
  )
)
```

### Profile Types
- `gene_profile()`: For annotated genes
- `align_profile()`: For read alignments  
- `axis_profile()`: For coordinate axes
- **Custom profiles**: Advanced users can create profiles for any genomic data type

### Profile Parameters
Each profile can define custom parameters that appear in the UI:
- **Select dropdowns**: For choosing between options
- **Numeric inputs**: For thresholds, sizes, etc.
- **Boolean toggles**: For show/hide features
- **Real-time updates**: Changes immediately affect visualization

## Data Requirements

Your data tables should follow these formats:

### Assembly Table
```
ASSEMBLY_ID	LIB_IDS
sample1	lib1,lib2,lib3
sample2	lib4,lib5,lib6
```

### Contig Table
```
contig	length	coverage	circular
contig_001	125000	5.2	FALSE
contig_002	87500	12.1	FALSE
contig_003	45000	3.8	TRUE
```

### Gene Table
```
gene	contig	start	end	strand	length	aa_length
gene_001	contig_001	100	1500	+	1401	467
gene_002	contig_001	1600	2800	-	1201	400
gene_003	contig_002	50	900	+	851	283
```

### Annotation Table (Optional)
```
gene	uniref	identity	coverage	evalue	prot_desc	tax
gene_001	UniRef100_A0A123	85.2	0.95	1e-150	Hypothetical protein	Bacteria
gene_002	UniRef100_B1B456	92.1	0.88	1e-200	DNA polymerase	Firmicutes
```

## Configuration Architecture

### Core Components
1. **Cache System**: Manages data loading and caching
2. **Data Module**: Handles lookup tables and data access
3. **Profile Manager**: Registers and manages visualization profiles
4. **Parameter System**: Handles UI controls and user preferences
5. **State Management**: Maintains application state across interactions

### File Structure
```
configs/your_config/
├── your_config_cfg.r          # Main configuration
│   ├── Data setup (lookup tables)
│   ├── Assembly/contig/genome registration
│   ├── View registration
│   └── Tab registration
└── your_config_main_view.r    # Profile definitions
    ├── Profile registration
    ├── Parameter definitions
    └── Data loading functions
```

### Customization Points
- **Data functions**: Modify how data is loaded and processed
- **Profile parameters**: Define user-controllable settings
- **Color schemes**: Customize visualization colors
- **Layout**: Adjust profile heights and arrangements
- **Interactions**: Add custom mouse/keyboard behaviors