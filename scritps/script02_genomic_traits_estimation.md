# Command lines and script for the computation of genomic traits for the chemostat experiment

### by Sara Beier & Angel Rain

The current file contains shell command lines for computing genomic traits as well from sequences output from the data2 as described in https://github.com/angelrainf/gesifus.cryo.dada2 and the genome data as given in Table S1 (https://github.com/sarabeier/genomic.traits) in the manuscript Beier et al, doi: https://doi.org/10.1101/2021.07.23.453341.


## Getting trait-specific tables 
For genome prediction by PICRUST, custom trait data was extracted from Table.S1 (https://github.com/sarabeier/genomic.traits). Individual datasets were formatted in R as follow and used for PICRUST2 prediction.

```R
library(tidyverse)

data.range<-function(n)
{
  # transform to range from 0-1000, round by two digits
  round((n-min(n,na.rm=T))/(max(n,na.rm=T)-min(n,na.rm=T))*1000,2)
}

gtraits <- tibble(read.table("TableS1.2021.tsv", header = T, sep = "\t", fill = T)) %>%
  select(IMG.Genome.ID, PICRUSt.ID,
         Genome.size, X.GC, RRN_IMG, RRN_rrnDB, X.HGT, CUB.F., Generation.time..Vieira.Silva., Generation.time..gRodon., Gene.duplication, Gene.richness, X.TF,
         Prophages) %>% 
  mutate(per.HGT.corr = ifelse(X.HGT < 65, X.HGT, NA)) %>%  # create column with %HGT values larger 65% 
  mutate(Prophages.corr = ifelse(Prophages < 20, Prophages, NA)) # create column with Prophages >2 
gtraits

ptraits <- data.frame(assembly=gtraits$PICRUSt.ID,
                      genome.size=data.range(gtraits$Genome.size),
                      HGT_perc=data.range(log(gtraits$per.HGT.corr+0.01)),
                      d=data.range(log(gtraits$Generation.time..Vieira.Silva.)),
                      F=data.range(gtraits$CUB.F.),
                      d.gRodon=data.range(log(gtraits$Generation.time..gRodon.)),
                      GC=data.range(gtraits$X.GC),
                      TF_perc=data.range(gtraits$X.TF),
                      mean.cog=data.range(log(gtraits$Gene.duplication)),
                      rich.cog=data.range(gtraits$Gene.richness),
                      prophages=data.range(gtraits$Prophages.corr))


var <-c('genomesize', 'HGT_perc', 
        'cub.d', 'cub.F',
        'TF_perc','GC',
        'cog.meancopynr','cog.rich',
        'vir', 'cub.dRg')

var.min <- c(min(gtraits[3],na.rm=T), min(gtraits[15],na.rm=T),
             min(gtraits[9],na.rm=T), min(gtraits[8],na.rm=T),
             min(gtraits[13],na.rm=T),min(gtraits[,4],na.rm=T),
             min(gtraits[11],na.rm=T), min(gtraits[12],na.rm=T),
             min(gtraits[16],na.rm=T), min(gtraits[,10],na.rm=T))

var.max <- c(max(gtraits[3],na.rm=T), max(gtraits[15],na.rm=T), 
             max(gtraits[9],na.rm=T), max(gtraits[8],na.rm=T),
             max(gtraits[13],na.rm=T),max(gtraits[,4],na.rm=T),
             max(gtraits[11],na.rm=T), max(gtraits[12],na.rm=T), 
             max(gtraits[16],na.rm=T), max(gtraits[,10],na.rm=T))

var.trans <- c('NA', 'log(x+0.01)', 
               'log(x)','NA',
               'NA','NA',
               'log(x)','NA',
               'NA','log(x)')

pic.traits <- data.frame(var, var.min, var.max,var.trans)
pic.traits
#Export transformation factors
write.table (pic.traits, 'data/picrust2.output/trait.data/pic.traits_trans.tab', sep='\t', row.names=FALSE)

#Export only traits associated to resilience and resistance traits as single files for picrust
#Resilience-related traits

write.table(na.omit(ptraits[c(1,6)]), 'data/picrust2.output/trait.data/p.d.gRodon.txt', row.names=FALSE, sep ='\t',quote=FALSE)
#Resistance-related traits
write.table(na.omit(ptraits[c(1,8)]), 'data/picrust2.output/trait.data/p.TF_perc.txt', row.names=FALSE, sep ='\t',quote=FALSE)
write.table(na.omit(ptraits[c(1,2)]), 'data/picrust2.output/trait.data/p.genome.size.txt', row.names=FALSE, sep ='\t',quote=FALSE)

```
Citation:  
* Beier,S., Werner,J., Bouvier,T., Mouquet,N., Violle, C. (2021). Trait-trait relationships and tradeoffs vary with genome size in prokaryotes. doi: https://doi.org/10.1101/2021.07.23.453341
* Douglas, G. M., Maffei, V. J., Zaneveld, J. R., Yurgel, S. N., Brown, J. R., Taylor, C. M., et al. (2020). PICRUSt2 for prediction of metagenome functions. Nature Biotechnology 38, 685–688. doi: 10.1038/s41587-020-0548-6.


## Trait-based estimations of gRodon, genome size, %TF

```bash
# Edit fasta.files sequences (Convert to multiline fastafiles)
sed -e 's/^/>/' -e 's/\ /\n/' dada-chem.seqs.nochim.fasta > dada-chem.seqs.nochim2.fasta

# PICRUST2 full pipeline
nohup picrust2_pipeline.py -s dada-chem.seqs.nochim2.fasta -i dada-chem.rcounts.ASV.2.tab --no_pathways --min_samples 0 -o picrust2_out_full -p 15 &

# Hidden state prediction (hsp) using custom trait tables
# generation time (gRodon)
hsp.py --observed_trait_table trait.data/p.d.gRodon.txt -t picrust_out_v1/out.tre -o picrust.traits/pic.dgR_predicted -n -p 10 #new

# Genome size
hsp.py --observed_trait_table trait.data/p.genome.size.txt -t picrust_out_v1/out.tre -o picrust.traits/pic.genomesize_predicted -n -p 10 #new
#%TF
hsp.py --observed_trait_table trait.data/p.TF_perc.txt -t picrust_out_v1/out.tre -o picrust.traits/pic.TF_perc_predicted -n -p 10 #new
```

# Prediction 16S copy number with custom tree

RRN prediction was performed by ausing the custom reference tree in PICRUST2 as described in https://github.com/picrust/picrust2/wiki/Sequence-placement and using the custom tree from Beier et al., 2021 ( https://github.com/sarabeier/genomic.traits/blob/main/scripts_phylogenetic.signal.md)

```bash
source activate picrust2
# a) Create custom references
# a1) rename muscle2.rrnDB.fasta
mv muscle2.rrnDB.fasta pro_ref.fasta
# a2) Rename rrndb_fasttree.format.tree
mv rrndb_fasttree.format.tree pro_ref.tre
# a3) Build a hmm profile for the aligned sequences
nohup hmmbuild pro_ref.hmm pro_ref.fasta &

# Check aligment sequences
raxml-ng --check --msa ../reference_tree_SB/pro_ref.fasta --model GTR+G --prefix T1

    #Evaluating the Model Parameters (https://github.com/Pbdas/epa-ng#setting-the-model-parameters)
    nohup raxml-ng --evaluate --msa ../reference_tree_SB/pro_ref.fasta --tree ../reference_tree_SB/pro_ref.tre --prefix info --model GTR+G+F --threads 10 &
    #Information required for picrust2 is store in info.raxml.bestModel> to be renamed > xxx.model
    mv info.raxml.bestModel /data/chemostats/input.dada2/picrust.with.updated.tree/custom_ref/custom_ref.model

    #custom_ref/
cd /data/chemostats/input.dada2/picrust.with.updated.tree/custom_ref/
cp pro_ref.fasta /data/chemostats/input.dada2/picrust.with.updated.tree/custom_ref/custom_ref.fasta
cp pro_ref.hmm /data/chemostats/input.dada2/picrust.with.updated.tree/custom_ref/custom_ref.hmm
cp pro_ref.tre /data/chemostats/input.dada2/picrust.with.updated.tree/custom_ref/custom_ref.tre
```

# 16S copy number prediction from custom tree
For prediction from custom reference tree PICRUST2 v2.4.2 was used since reference tree application was implemented after version 2.4.0. For our trait prediction RRN mean values (16S_SB.txt) was used as computed by Beier et al., 2021 using the extrapolated RRN values from the rrnDB database (rrnDB, Stoddard et al. 2015)

```bash
# Run picrust pipeline with custom tree
## 16S copy number trait from Beier et al., 2021 (see above)
16S_SB.txt
## 16S copy number
cd /data/chemostats/input.dada2/picrust2_out_custom_treev2

nohup picrust2_pipeline.py -s ../dada-chem.seqs.nochim2.fasta -i ../dada-chem.rcounts.ASV.2.tab --min_samples 0 -o picrust2_out_custom_treev2 -p 10 --ref_dir ../picrust.with.updated.tree/custom_ref/ --marker_gene_table ../picrust.with.updated.tree/trait_tables/16S_SB.txt --no_pathways --skip_minpath --edge_exponent 0 &

# Unzip file
gunzip marker_predicted_and_nsti.tsv.gz
# Renaming it for downstream analysis in R
mv marker_predicted_and_nsti.tsv pic.16S_predicted_custom_tree.txt 
```
Citation:  
* Beier,S., Werner,J., Bouvier,T., Mouquet,N., Violle, C. (2022). Trait-trait relationships and tradeoffs vary with genome size in prokaryotes. Front. Microbiol. doi: https://doi.org/10.3389/fmicb.2022.985216
* Stoddard, S. F., Smith, B. J., Hein, R., Roller, B. R. K., and Schmidt, T. M. (2015). rrnDB: improved tools for interpreting rRNA gene abundance in bacteria and archaea and a new foundation for future development. Nucleic Acids Res. 43, D593–D598. doi:10.1093/nar/gku1201.


## Re-transform trait-specific tables 
For genome prediction by PICRUST2, trait data were re-transformed  back to the raw units in R as follow:

```R
# PICRUST data re-transformation
#Loading convertion factors
pic.traits <- read.table("../data/picrust2/trait.data/pic.traits_trans.tab", header=T, sep ='\t')

#load generation time predicted from PICRUSt2 default settings
pic.d.gRodon.default <- read.table("../data/picrust2/trait.predicted/pic.chemo10.dgR_predicted.txt", header=T)  

#load %TF predicted from PICRUSt2 default settings
pic.TFr.default <- read.table("../data/picrust2/trait.predicted/pic.chemo10.TF_perc_predicted.txt", header=T)

#load genome size predicted from PICRUSt2 default settings
pic.gs.default <- read.table("../data/picrust2/trait.predicted/pic.chemo10.genomesize_predicted.txt", header=T)

#Re-transform predicted genomic traits

## Generation time
n <- pic.d.gRodon.default$d.gRodon
pic.d.gRodon.default$d.gRodon <- exp(1)^ (n / 1000 * (log(pic.traits$var.max[10])-log(pic.traits$var.min[10])) + log(pic.traits$var.min[10])) #log(x) transformation

## Transcription factors
n <- pic.TFr.default$TF_perc
pic.TFr.default$TF_perc <- (n / 1000 * (pic.traits$var.max[5]-pic.traits$var.min[5])) + pic.traits$var.min[5] 

## Genome size
n <- pic.gs.default$genome.size
pic.gs.default$genome.size <- ((n / 1000 * (pic.traits$var.max[1]-pic.traits$var.min[1])) + pic.traits$var.min[1])/1000000 #divide by MiO to get results in Mbp

Exporting re-transformed traits
write.table(pic.d.gRodon.default, '../data/picrust2/trait.predicted/pic.d.gRodon.scaled.txt', row.names=FALSE, sep ='\t',quote=FALSE)
write.table(pic.TFr.default, '../data/picrust2/trait.predicted/pic.TF_perc.scaled.txt', row.names=FALSE, sep ='\t',quote=FALSE)
write.table(pic.gs.default, '../data/picrust2/trait.predicted/pic.genome.size.scaled.txt', row.names=FALSE, sep ='\t',quote=FALSE)
```

