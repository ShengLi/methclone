methclone: Detect the dynamic evolution of clonal epialleles in DNA methylation sequencing data.
version: 0.1

==Installation==
cd /path/to/methclone/src/utils/BamTools/
mkdir -p lib
make
cd /path/to/methclone/src/utils/gzstream/
make 
cd /path/to/methclone/src/
mkdir -p ../bin
make

==Usage==
./methclone stage1.bam stage2.bam output.txt.gz sampleID

==Example==
cd /path/to/methclone/
./bin/methclone example/chr22-stage1.bam example/chr22-stage2.bam example/chr22.output.txt.gz chr22-1-2
