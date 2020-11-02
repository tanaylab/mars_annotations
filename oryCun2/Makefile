# UCSC_ROOT := http://hgdownload.soe.ucsc.edu/gbdb/oryCun2/ncbiRefSeq/
UCSC_ROOT := http://hgdownload.soe.ucsc.edu/goldenPath/oryCun2/database/

.PHONY: all clean

all: gene_intervals_oryCun2.txt.gz

clean:
	\rm -rf gene_intervals_oryCun2.txt
	\rm -rf ncbiRefSeq.txt.gz
	\rm -rf ncbiRefSeqOther.bb

gene_intervals_oryCun2.txt.gz: ncbiRefSeq.txt.gz
	./build_intervals.R  $<  $@

ncbiRefSeq.txt.gz:
	wget $(UCSC_ROOT)/$@
