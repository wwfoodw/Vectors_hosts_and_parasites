mkdir ~/may_21
cd ~/may_21

module load Parallel/20210422

gunzip_func() {
for f in $1
do
gunzip "$f"
done
}
# export function
export -f gunzip_func

## run function in parallel
# $1 filenames to use in parallel,

parallel -j 8 --line-buffer gunzip_func ::: *.gz

###################
#                 #
#### MOSQUITO #####
#                 #
###################

module load cutadapt/3.3-gimkl-2020a-Python-3.8.2

fastq1=($(find *_R1.fastq))
fastq2=($(find *_R2.fastq))
long=($(echo ${#fastq1[@]}))
for i in `seq 0 $((long-1))`
do
# i5 is g
# i7 is G
cutadapt --cores=0 --action=trim --discard-untrimmed -g ^CATCCAATCCATAATAAAGCAT -G ^AAAAATACCCTTCTATCCAAATCT -o ant_mal-${fastq1[i]%%.*}_cutadapt.fastq -p ant_mal-${fastq2[i]%%.*}_cutadapt.fastq ${fastq1[i]} ${fastq2[i]}
cutadapt --cores=0 --action=trim --discard-untrimmed -g ^GGGCAATCCTGAGCCAA -G ^CCATTGAGTCTCTGCACCTATC -o trnl-${fastq1[i]%%.*}_cutadapt.fastq -p trnl-${fastq2[i]%%.*}_cutadapt.fastq ${fastq1[i]} ${fastq2[i]}
cutadapt --cores=0 --action=trim --discard-untrimmed -g ^TNTTYTCMACYAACCACAAAGA -G ^CARAAGCTYATGTTRTTYATDCG -o rep-${fastq1[i]%%.*}_cutadapt.fastq -p rep-${fastq2[i]%%.*}_cutadapt.fastq ${fastq1[i]} ${fastq2[i]}
cutadapt --cores=0 --action=trim --discard-untrimmed -g ^GGWACWGGWTGAACWGTWTAYCCYCC -G ^TANACYTCNGGRTGNCCRAARAAYCA -o coi-${fastq1[i]%%.*}_cutadapt.fastq -p coi-${fastq2[i]%%.*}_cutadapt.fastq ${fastq1[i]} ${fastq2[i]}
done

find . -type f -empty -print -delete # remove empty files

mkdir ~/may_21/mosquito

cp ant_mal-* mosquito/
cp trnl-* mosquito/
cp rep-* mosquito/
cp coi-* mosquito/

cd ~/may_21/mosquito

#############  #############  #############  #############  #############  #############  #############  #############
#
# # subsample these down to 100,000
#
#############  #############  #############  #############  #############  #############  #############  #############

cd ~/may_21/mosquito

# make a list of file sizes
rm count_fastq.txt
rm large_fastq.txt
rm small_fastq.txt
rm remove_fastq.txt
for i in *.fastq
do
    wc "$i" -l >> count_fastq.txt
done

# make a list of files with > 400,000 lines (100,000 reads)
awk '{if($1>400000)print$2}' count_fastq.txt > large_fastq.txt
awk '{if($1<400001)print$2}' count_fastq.txt > small_fastq.txt

module load SeqKit/2.2.0
module load seqtk/1.3-gimkl-2018b

# subsample these down to 100,000 
subsample_func() {
for f_file in $1
do
seqtk sample -s123 "$f_file" 100000 > "$f_file.subsample.fastq" # note: -s is the SEED, if doing paired ends #20000
done
}
export -f subsample_func # export function

parallel -j 8 --line-buffer subsample_func :::: large_fastq.txt ## run function in parallel

rename_func() {
for f_file in $1
do
cp "$f_file" "$f_file.subsample.fastq"
done
}
export -f rename_func # export function

parallel -j 8 --line-buffer rename_func :::: small_fastq.txt ## run function in parallel

mkdir ant_mal
mkdir trnl
mkdir rep
mkdir coi

mv ant_mal*.subsample.fastq ant_mal/
mv trnl*.subsample.fastq trnl/
mv rep*.subsample.fastq rep/
mv coi*.subsample.fastq coi/

#############  #############  #############  #############  #############  #############  #############  #############
#
# Dada2 coi
#
#############  #############  #############  #############  #############  #############  #############  #############

cd ~/may_21/mosquito

module load  GSL/2.6-GCCcore-9.2.0
module load  R/4.1.0-gimkl-2020a

R

## Load the necessary libraries**
library("dada2")
library("ggplot2")
library("Biostrings")

## Set up directory

rm(list=ls())
set.seed(123)

miseq_path <- file.path("coi/")
  
# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(miseq_path, pattern="_L007_R1_cutadapt.fastq.subsample.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_L007_R2_cutadapt.fastq.subsample.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- paste(sapply(strsplit(fnFs, "_"), `[`, 4), sep="_")
sampleNames
length(sampleNames)
sampleNames[duplicated(sampleNames)] # check for duplicates
sampleNames <- paste(sapply(strsplit(fnFs, "_"), `[`, 4), sep="_")
sampleNames
length(sampleNames)
sampleNames[duplicated(sampleNames)] # check for duplicates

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)

filt_path <- file.path(miseq_path, "filtered_c") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))


# These should match:
length(fnFs)
length(fnRs)
length(filtFs)
length(filtRs)
length(sampleNames)

# Filter the forward and reverse reads:

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  maxN=0,
  maxEE=c(2,2),
  truncQ=10,
  rm.phix=TRUE,
  compress=FALSE,
  multithread=FALSE
)
      
head(out)

ob <- as.data.frame(out)
new_out <- ob[ob$reads.out == 0, ]

# Dereplicate

filtFs <- file.path(filt_path, (sort(list.files(filt_path, pattern="_F_filt.fastq.gz"))))
filtRs <- file.path(filt_path, (sort(list.files(filt_path, pattern="_R_filt.fastq.gz"))))

derepFs <- derepFastq(filtFs, verbose=FALSE)
derepRs <- derepFastq(filtRs, verbose=FALSE)

# Name the derep-class objects by the sample names
x <- sort(list.files(filt_path, pattern="_F_filt.fastq.gz"))
sampleNames <- paste(sapply(strsplit(x, "_"), `[`, 1), sapply(strsplit(x, "_"), `[`, 3), sep="_")
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)


dadaFs <- dada(derepFs, err=errF, pool = TRUE, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, pool = TRUE, multithread=TRUE)

# Inspect the dada-class object
dadaFs[[1]]
dadaRs[[1]]

# Construct sequence table
cat_mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, justConcatenate=TRUE) # DON'T MERGE - Concatenate
cat_seqtabAll <- makeSequenceTable(cat_mergers[!grepl("Mock", names(cat_mergers))])
table(nchar(getSequences(cat_seqtabAll)))

cat_seqtabAll <- cat_seqtabAll[, colSums(cat_seqtabAll) > 100]
cat_seqtabAll <- cat_seqtabAll[ rowSums(cat_seqtabAll) > 100, ]

coi_cat_sequences <- colnames(cat_seqtabAll)
colnames(cat_seqtabAll) <- paste("Hap_", 1:length(colnames(cat_seqtabAll)), sep ="")
names(coi_cat_sequences) <- colnames(cat_seqtabAll)

ls()
rm("dadaFs","dadaRs","derepFs","derepRs",
    "errF","errR","filt_path","filtFs",
    "filtRs","fnFs","fnRs","mergers",
    "miseq_path","new_out","ob","out",
    "sampleNames","x","seqtabAll","cat_mergers","coi_sequences")

writeXStringSet(DNAStringSet(coi_sequences), file = "ps_coi.fasta")
writeXStringSet(DNAStringSet(coi_cat_sequences), file = "ps_cat_coi.fasta")
save.image(file = "dada2_coi.RData")

#############  #############  #############  #############  #############  #############  #############  #############
#
# Dada2 trnl
#
#############  #############  #############  #############  #############  #############  #############  #############

cd ~/may_21/mosquito

module load  GSL/2.6-GCCcore-9.2.0
module load  R/4.1.0-gimkl-2020a

R

## Load the necessary libraries**
library("dada2")
library("ggplot2")
library("Biostrings")

## Set up directory

rm(list=ls())
set.seed(123)

miseq_path <- file.path("trnl/")
  
# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(miseq_path, pattern="_L007_R1_cutadapt.fastq.subsample.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_L007_R2_cutadapt.fastq.subsample.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- paste(sapply(strsplit(fnFs, "_L007"), `[`, 1),
                            sep="_")

length(sampleNames)
sampleNames[duplicated(sampleNames)] # check for duplicates

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)

filt_path <- file.path(miseq_path, "filtered_c") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))


# These should match:
length(fnFs)
length(fnRs)
length(filtFs)
length(filtRs)
length(sampleNames)

# Filter the forward and reverse reads:

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  maxN=0,
  maxEE=c(2,2),
  truncQ=10,
  rm.phix=TRUE,
  compress=FALSE,
  multithread=FALSE
)
      
head(out)

ob <- as.data.frame(out)
new_out <- ob[ob$reads.out == 0, ]

# Dereplicate

filtFs <- file.path(filt_path, (sort(list.files(filt_path, pattern="_F_filt.fastq.gz"))))
filtRs <- file.path(filt_path, (sort(list.files(filt_path, pattern="_R_filt.fastq.gz"))))

derepFs <- derepFastq(filtFs, verbose=FALSE)
derepRs <- derepFastq(filtRs, verbose=FALSE)

# Name the derep-class objects by the sample names
x <- sort(list.files(filt_path, pattern="_F_filt.fastq.gz"))
sampleNames <- paste(sapply(strsplit(x, "_"), `[`, 1), sapply(strsplit(x, "_"), `[`, 3), sep="_")
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)


dadaFs <- dada(derepFs, err=errF, pool = TRUE, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, pool = TRUE, multithread=TRUE)

# Inspect the dada-class object
dadaFs[[1]]
dadaRs[[1]]


# Construct sequence table
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs) 
seqtabAll <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
table(nchar(getSequences(seqtabAll)))

apply(seqtabAll,1,function(x){sum(x>0)})

trnl_sequences <- colnames(seqtabAll)
colnames(seqtabAll) <- paste("Hap_", 1:length(colnames(seqtabAll)), sep ="")
names(trnl_sequences) <- colnames(seqtabAll)

# Construct sequence table
cat_mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, justConcatenate=TRUE) # DON'T MERGE - Concatenate
cat_seqtabAll <- makeSequenceTable(cat_mergers[!grepl("Mock", names(cat_mergers))])
table(nchar(getSequences(cat_seqtabAll)))

apply(cat_seqtabAll,1,function(x){sum(x>0)})

trnl_cat_sequences <- colnames(cat_seqtabAll)
colnames(cat_seqtabAll) <- paste("Hap_", 1:length(colnames(cat_seqtabAll)), sep ="")
names(trnl_cat_sequences) <- colnames(cat_seqtabAll)

ls()
rm("dadaFs","dadaRs","derepFs","derepRs",
    "errF","errR","filt_path","filtFs",
    "filtRs","fnFs","fnRs","mergers",
    "miseq_path","new_out","ob","out",
    "sampleNames","x")

writeXStringSet(DNAStringSet(trnl_sequences), file = "ps_trnl.fasta")
writeXStringSet(DNAStringSet(trnl_cat_sequences), file = "ps_cat_trnl.fasta")
save.image(file = "dada2_trnl.RData")


#############  #############  #############  #############  #############  #############  #############  #############
#
# Dada2 rep
#
#############  #############  #############  #############  #############  #############  #############  #############

cd ~/may_21/mosquito

module load  GSL/2.6-GCCcore-9.2.0
module load  R/4.1.0-gimkl-2020a

R

## Load the necessary libraries**
library("dada2")
library("ggplot2")
library("Biostrings")

## Set up directory

rm(list=ls())
set.seed(123)

miseq_path <- file.path("rep/")
  
# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(miseq_path, pattern="_L007_R1_cutadapt.fastq.subsample.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_L007_R2_cutadapt.fastq.subsample.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- paste(sapply(strsplit(fnFs, "_L007"), `[`, 1),
                            sep="_")

length(sampleNames)
sampleNames[duplicated(sampleNames)] # check for duplicates

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)

filt_path <- file.path(miseq_path, "filtered_c") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))


# These should match:
length(fnFs)
length(fnRs)
length(filtFs)
length(filtRs)
length(sampleNames)

# Filter the forward and reverse reads:

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  maxN=0,
  maxEE=c(2,2),
  truncQ=10,
  rm.phix=TRUE,
  compress=FALSE,
  multithread=FALSE
)
      
head(out)

ob <- as.data.frame(out)
new_out <- ob[ob$reads.out == 0, ]

# Dereplicate

filtFs <- file.path(filt_path, (sort(list.files(filt_path, pattern="_F_filt.fastq.gz"))))
filtRs <- file.path(filt_path, (sort(list.files(filt_path, pattern="_R_filt.fastq.gz"))))

derepFs <- derepFastq(filtFs, verbose=FALSE)
derepRs <- derepFastq(filtRs, verbose=FALSE)

# Name the derep-class objects by the sample names
x <- sort(list.files(filt_path, pattern="_F_filt.fastq.gz"))
sampleNames <- paste(sapply(strsplit(x, "_"), `[`, 1), sapply(strsplit(x, "_"), `[`, 3), sep="_")
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)


dadaFs <- dada(derepFs, err=errF, pool = TRUE, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, pool = TRUE, multithread=TRUE)

# Inspect the dada-class object
dadaFs[[1]]
dadaRs[[1]]


# Construct sequence table
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs) 
seqtabAll <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
table(nchar(getSequences(seqtabAll)))

apply(seqtabAll,1,function(x){sum(x>0)})

rep_sequences <- colnames(seqtabAll)
colnames(seqtabAll) <- paste("Hap_", 1:length(colnames(seqtabAll)), sep ="")
names(rep_sequences) <- colnames(seqtabAll)

# Construct sequence table
cat_mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, justConcatenate=TRUE) # DON'T MERGE - Concatenate
cat_seqtabAll <- makeSequenceTable(cat_mergers[!grepl("Mock", names(cat_mergers))])
table(nchar(getSequences(cat_seqtabAll)))

apply(cat_seqtabAll,1,function(x){sum(x>0)})

rep_cat_sequences <- colnames(cat_seqtabAll)
colnames(cat_seqtabAll) <- paste("Hap_", 1:length(colnames(cat_seqtabAll)), sep ="")
names(rep_cat_sequences) <- colnames(cat_seqtabAll)

ls()
rm("dadaFs","dadaRs","derepFs","derepRs",
    "errF","errR","filt_path","filtFs",
    "filtRs","fnFs","fnRs","mergers",
    "miseq_path","new_out","ob","out",
    "sampleNames","x")

writeXStringSet(DNAStringSet(rep_sequences), file = "ps_rep.fasta")
writeXStringSet(DNAStringSet(rep_cat_sequences), file = "ps_cat_rep.fasta")
save.image(file = "dada2_rep.RData")


#############  #############  #############  #############  #############  #############  #############  #############
#
# Dada2 ant_mal
#
#############  #############  #############  #############  #############  #############  #############  #############

cd ~/may_21/mosquito

module load  GSL/2.6-GCCcore-9.2.0
module load  R/4.1.0-gimkl-2020a

R

## Load the necessary libraries**
library("dada2")
library("ggplot2")
library("Biostrings")

## Set up directory

rm(list=ls())
set.seed(123)

miseq_path <- file.path("ant_mal/")
  
# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(miseq_path, pattern="_L007_R1_cutadapt.fastq.subsample.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_L007_R2_cutadapt.fastq.subsample.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- paste(sapply(strsplit(fnFs, "_L007"), `[`, 1),
                            sep="_")

length(sampleNames)
sampleNames[duplicated(sampleNames)] # check for duplicates

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)

filt_path <- file.path(miseq_path, "filtered_c") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))


# These should match:
length(fnFs)
length(fnRs)
length(filtFs)
length(filtRs)
length(sampleNames)

# Filter the forward and reverse reads:

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  maxN=0,
  maxEE=c(2,2),
  truncQ=10,
  rm.phix=TRUE,
  compress=FALSE,
  multithread=FALSE
)
      
head(out)

ob <- as.data.frame(out)
new_out <- ob[ob$reads.out == 0, ]

# Dereplicate

filtFs <- file.path(filt_path, (sort(list.files(filt_path, pattern="_F_filt.fastq.gz"))))
filtRs <- file.path(filt_path, (sort(list.files(filt_path, pattern="_R_filt.fastq.gz"))))

derepFs <- derepFastq(filtFs, verbose=FALSE)
derepRs <- derepFastq(filtRs, verbose=FALSE)

# Name the derep-class objects by the sample names
x <- sort(list.files(filt_path, pattern="_F_filt.fastq.gz"))
sampleNames <- paste(sapply(strsplit(x, "_"), `[`, 1), sapply(strsplit(x, "_"), `[`, 3), sep="_")
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)


dadaFs <- dada(derepFs, err=errF, pool = TRUE, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, pool = TRUE, multithread=TRUE)

# Inspect the dada-class object
dadaFs[[1]]
dadaRs[[1]]


# Construct sequence table
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs) 
seqtabAll <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
table(nchar(getSequences(seqtabAll)))

apply(seqtabAll,1,function(x){sum(x>0)})

ant_mal_sequences <- colnames(seqtabAll)
colnames(seqtabAll) <- paste("Hap_", 1:length(colnames(seqtabAll)), sep ="")
names(ant_mal_sequences) <- colnames(seqtabAll)

# Construct sequence table
cat_mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, justConcatenate=TRUE) # DON'T MERGE - Concatenate
cat_seqtabAll <- makeSequenceTable(cat_mergers[!grepl("Mock", names(cat_mergers))])
table(nchar(getSequences(cat_seqtabAll)))

apply(cat_seqtabAll,1,function(x){sum(x>0)})

ant_mal_cat_sequences <- colnames(cat_seqtabAll)
colnames(cat_seqtabAll) <- paste("Hap_", 1:length(colnames(cat_seqtabAll)), sep ="")
names(ant_mal_cat_sequences) <- colnames(cat_seqtabAll)

ls()
rm("dadaFs","dadaRs","derepFs","derepRs",
    "errF","errR","filt_path","filtFs",
    "filtRs","fnFs","fnRs","mergers",
    "miseq_path","new_out","ob","out",
    "sampleNames","x")

writeXStringSet(DNAStringSet(ant_mal_sequences), file = "ps_ant_mal.fasta")
writeXStringSet(DNAStringSet(ant_mal_cat_sequences), file = "ps_cat_ant_mal.fasta")
save.image(file = "dada2_ant_mal.RData")


###
#
# process fasta files
#
###

cd ~/may_21/mosquito

awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq=seq $0} END {if (seq) print seq}' \
  ps_cat_trnl.fasta > trnl_single_line.fasta
awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq=seq $0} END {if (seq) print seq}' \
  ps_cat_coi.fasta > coi_single_line.fasta
awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq=seq $0} END {if (seq) print seq}' \
  ps_cat_rep.fasta > rep_single_line.fasta
awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq=seq $0} END {if (seq) print seq}' \
  ps_cat_ant_mal.fasta > ant_mal_single_line.fasta


sed 's/^.*NNNNNNNN//' coi_single_line.fasta > coi_right_single_line.fasta
sed 's/^.*NNNNNNNN//' ant_mal_single_line.fasta > ant_right_mal_single_line.fasta

sed 's/NNNNNNNN.*$//' coi_single_line.fasta > coi_left_single_line.fasta
sed 's/NNNNNNNN.*$//' ant_mal_single_line.fasta > ant_left_mal_single_line.fasta

###
###
###


