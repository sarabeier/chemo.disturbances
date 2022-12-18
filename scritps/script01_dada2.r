#please log into the bio40 server: ssh -X sbeier@10.11.20.40

#the path has slightly changed there is no folder /bio49 anymore but instead only bio

#the best is to switch to phy-1 (ssh phy-1) and start there R by enteringe module load R

#create a new folder for output data /bio/projects/2019/gesifus.chemostats/analyses/dada.out4/

#chem.dada4.R
#analyzes sequence data from Delivery.1 and Delivery.2 in folder 'rawdata/metabercoding/PrimerClipped_joined' condsidered
#output folder: 'analyses/dada.out4/'
#created by Sara Beier

print('run chem.dada2.R')

rm(list=ls())
library(microbiome)
library(dada2)
library(phyloseq)
library(ggplot2)
library(gridExtra)
library(DECIPHER)
library(plyr)
library(phangorn)
library(vegan)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
print('Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq')
path<-("/bio/projects/2019/gesifus.chemostats/rawdata/metabarcoding/PrimerClipped_joined")
list.files(path)

fnFs <- sort(list.files(path, pattern=glob2rx("*-C10-*_R1.fastq"), full.names = TRUE)) #SARA: only C10
fnRs <- sort(list.files(path, pattern=glob2rx("*-C10-*_R2.fastq"), full.names = TRUE)) #SARA: only C10

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq and create dataframe with metadata
print('Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq and create dataframe with metadata')
sample.ID <- substring(sapply(strsplit(basename(fnFs), "_"), `[`, 1),23,30)
T <- substring(sample.ID,5,5)
Chem.ID <-substring(sample.ID,7,8)
dat <- data.frame(sample.ID, Chem.ID, T)
head(dat, 20)


# Place filtered files in filtered/ subdirectory
print('Place filtered files in filtered/ subdirectory')
filtFs <- file.path(path, "filtered.chem", paste0(sample.ID, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered.chem", paste0(sample.ID, "_R_filt.fastq.gz"))
## the truncation length (truncLen) is based on what you see with `plotQualityProfile`
#all other arguments are on their default levels
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,180),orient.fwd ="TACG",maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)
#option orient.fwd ="TACG" decreases number of passed reads!!

#Learn the error rates of each base
print('Learn the error rates of each base')
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=FALSE) #multithread = FALSE to avoid error


#Remove replications
print('Remove replications')
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
print('Name the derep-class objects by the sample names')
names(derepFs) <- sample.ID
names(derepRs) <- sample.ID

#Sample inference
print('Sample inference')
dadaFs <- dada(derepFs, err=errF, multithread=FALSE, pool = TRUE) #multithread = FALSE to avoid error
dadaRs <- dada(derepRs, err=errR, multithread=FALSE, pool = TRUE) #multithread = FALSE to avoid error
dadaFs[[1]]
dadaRs[[1]]
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
print('Inspect the merger data.frame from the first sample')
head(mergers[[1]])

#Construct sequence table
print('Construct sequence table')
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
print('Inspect distribution of sequence lengths')
table(nchar(getSequences(seqtab)))

#remove chimeras
print('remove chimeras')
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE) #multithread = FALSE to avoid error
dim(seqtab.nochim)

print(paste (sum(seqtab.nochim)/sum(seqtab),': fraction of recovered reads after removal of chimeras'))
#We lost a lot of sequences. But when we account for the abundances of
#those variants we see they account for only about 15% of the merged sequence reads

#Track the number of reads that made it through so far
print('Track the number of reads that made it through so far')
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.ID
head(track)
write.table(track,"/chemostats/analyses/dada2.out/dada-chem.Number_of_reads_per_sample.tab", sep = "\t")

#Assign taxonomy. You have to download the `silva_nr_v128_train_set.fa.gz` file
#and put it inside the sae folder of your samples. See DADA2 documentation to download it
#remember to only use multithread=TRUE on Mac, not with Windows based computers
print('Assign taxonomy')
taxa.si <- assignTaxonomy(seqtab.nochim, "/db/dada2/silva_nr_v128_train_set.fa.gz",  multithread=F) #set multithread=F instead if multithread=TRUE to avoid error
taxa.gt <- assignTaxonomy(seqtab.nochim, "",  multithread=F) #set multithread=F instead if multithread=TRUE to avoid error

#Inspect taxonomic assignments
print('Inspect taxonomic assignments')
taxa.print <- taxa.gt # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
taxa.print.names<-as.data.frame(taxa.print)

#############################
#simple data.frame construction from the information encoded in the filenames
#Sample names are the 16S sequences -->ps0
samples.out <- rownames(seqtab.nochim)
#subject <- sapply(strsplit(samples.out, "D"), `[`, 1) #usuful to see your samples and set the metadata (in my case Community, DOM and Treatment)
samdf <- data.frame(sample.ID=sample.ID, Chem.ID=Chem.ID,T=T)
rownames(samdf) <- samples.out
##phylosq object from dada2 output
ps0 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
sample_data(samdf),
tax_table(taxa.gt))
ps0 <- prune_samples(sample_names(ps0) != "Mock", ps0) # Remove mock sample (if you have one)
ps0
#############################


#alignment of sequences (if you exported, can be done in, e.g., the Geneious program)
print('alignment of sequences and extract fasta file')
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree. This could be a problem because
#the default seqs names are the actual 16S sequence. Therefore, you might need to change these names
names(seqs)<- paste0("SV_", seq(ntaxa(ps0)), "_", taxa.print.names$Order)
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
seqs.tab<-as.data.frame(seqs)
write.table(seqs, "/chemostats/analyses/dada2.out/dada-chem.seqs.nochim.fasta",quote=FALSE, col.names=FALSE)


#simple data.frame construction from the information encoded in the filenames
#here I changed the sample names so I don't have the 16S sequences as names -->ps
seqtab.nochim2<-seqtab.nochim
colnames(seqtab.nochim2)<-rownames(seqs.tab)
rownames(taxa.gt)<-rownames(seqs.tab)
rownames(taxa.si)<-rownames(seqs.tab)


#I create another phyloseq object with the data with the new sample names.
#same procedure as before
print('phyloseq object with the data with the new sample names')
samples.out2 <- rownames(seqtab.nochim)
samdf <- data.frame(sample.ID=sample.ID, Chem.ID=Chem.ID, T=T)
rownames(samdf) <- samples.out2
#phylosq object from dada2 output
#with `seqtab.nochim2`, with the changed names (not the 16S sequences)
ps <- phyloseq(otu_table(seqtab.nochim2, taxa_are_rows=FALSE),
sample_data(samdf),
tax_table(taxa.gt))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps
saveRDS(ps, "/bio/projects/2019/gesifus.chemostats/analyses/dada.out4/chem.ps.rds")

#extract count table (raw count data)
print('extract count table (raw count data)')
OTU1 = as(otu_table(ps), "matrix") #temporary file
# transpose if necessary
if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf.ps = as.data.frame(OTU1)
OTUdf.ps <- OTUdf.ps[with(OTUdf.ps, order(row.names(OTUdf.ps))), ]


#extract relative abundance table
print('extract relative abundance table')
ps.rel = transform_sample_counts(ps, function(x) x/sum(x))
rOTU1 = as(otu_table(ps.rel), "matrix") #temporary file
# transpose if necessary
if(taxa_are_rows(ps.rel)){rOTU1 <- t(rOTU1)}
# Coerce to data.frame
rOTUdf.ps = as.data.frame(rOTU1)
rOTUdf.ps <- rOTUdf.ps[with(rOTUdf.ps, order(row.names(rOTUdf.ps))), ]
trOTUdf.ps <- as.data.frame.matrix(t(rOTUdf.ps))

#rarefy data
print('rarefy data')
OTUdf.rar.ps <-rrarefy(OTUdf.ps, min (rowSums(OTUdf.ps)))
rOTUdf.rar.ps <-prop.table(OTUdf.rar.ps, margin = 1)
trOTUdf.rar.ps <- as.data.frame.matrix(t(rOTUdf.rar.ps))

#write count-tables
print('write count-tables')

write.table(data.frame("ASV"=rownames(trOTUdf.ps),trOTUdf.ps),"/bio/projects/2019/gesifus.chemostats/analyses/dada.out4/dada-chem.rcounts.ASV.2.tab", row.names=FALSE,quote=FALSE, sep='\t') #I expert also the transposed data in this case, beause I need it like this for the picrust

##construction of a phylogenetic tree
print('construction of a phylogenetic tree')
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR2 <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

write.tree(fitGTR2$tree, file = "/chemostats/analyses/dada2.out/dada-chem.GTR2") #export tree
