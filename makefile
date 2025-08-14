create_minimal:
	mkdir -p configs/minimal
	Rscript utils/make_minimal.r	

# commands used to generate paf files for BAH and BAA, run within pg20
# cat `find /makeshift-mnt/output/minimap/v1.08/map-hifi/lib/BAH -name align.paf` | grep -P '\t(ctg1212|ctg3177|ctg520|ctg6585|ctg4261|ctg4393|ctg2302|ctg5193|ctg6760|ctg5008)\t' > /tmp/BAH.paf
# cat `find /makeshift-mnt/output/minimap/v1.08/map-hifi/lib/BAA -name align.paf` | grep -P '\t(ctg965|ctg2477|ctg1325|ctg1360|ctg404|ctg923|ctg1331|ctg986|ctg1956|ctg1384)\t' > /tmp/BAA.paf

AID ?= BAH
AIDS ?= BAA BAH
INPUT_PAF ?= examples/paf/$(AID).paf
OUTPUT_ALN ?= examples/aln/$(AID).aln
aln_construct:
	alntools construct -ifn_paf $(INPUT_PAF) -ofn $(OUTPUT_ALN)
aln_construct_all:
	mkdir -p examples/aln
	$(foreach aid, $(AIDS), $(MAKE) aln_construct AID=$(aid);)