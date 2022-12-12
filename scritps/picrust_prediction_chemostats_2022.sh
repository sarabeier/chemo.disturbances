#!/bin/sh
#  picrust_prediction_chemostats_2022.sh
#  Created by Angel on 19.11.22.
#  
# Accessing to the cloud
ssh -i /Users/angel/.ssh/cloud.key ubuntu@193.196.20.91
#Enter passphrase for key '/Users/angel/.ssh/cloud.key': ---> "Bonjour"

#Go to the directory
cd /data/chemostats/input.dada2_v2
# edit fasta.files sequences (Convert to multiline fastafiles)
sed -e 's/^/>/' -e 's/\ /\n/' dada-chem.seqs.nochim.fasta > dada-chem.seqs.nochim2.fasta
## Trait-based estimations of gRodon, genome size, %TF
# PICRUST2 full pipeline (version 2.3.0_b)
source activate picrust2

nohup picrust2_pipeline.py -s dada-chem.seqs.nochim2.fasta -i dada-chem.rcounts.ASV.2.tab --no_pathways --min_samples 0 -o picrust2_out_full -p 20 &

# Hidden state prediction (hsp) using custom trait tables
# generation time (gRodon)
hsp.py --observed_trait_table trait.files.new/p.d.gRodon.txt -t picrust2_out_full/out.tre -o picrust.traits/pic.dgR_predicted -n -p 10 #new

# Genome size
hsp.py --observed_trait_table trait.files.new/p.genome.size.txt -t picrust2_out_full/out.tre -o picrust.traits/pic.genomesize_predicted -n -p 10 #new
#%TF
hsp.py --observed_trait_table trait.files.new/p.TF_perc.txt -t picrust2_out_full/out.tre -o picrust.traits/pic.TF_perc_predicted -n -p 10 #new


# Prediction 16S copy number with custom tree

RRN prediction was performed by ausing the custom reference tree in PICRUST2 as described in https://github.com/picrust/picrust2/wiki/Sequence-placement and using the custom tree from Beier et al., 2021 ( https://github.com/sarabeier/genomic.traits/blob/main/scripts_phylogenetic.signal.md)


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

# 16S copy number prediction from custom tree
For prediction from custom reference tree PICRUST2 v2.4.2 was used since reference tree application was implemented after version 2.4.0. For our trait prediction RRN mean values (16S_SB.txt) was used as computed by Beier et al., 2021 using the extrapolated RRN values from the rrnDB database (rrnDB, Stoddard et al. 2015)

# Run picrust pipeline with custom tree
## 16S copy number trait from Beier et al., 2021 (see above)
16S_SB.txt
## 16S copy number
cd /data/chemostats/input.dada2/picrust2_out_custom_treev2
conda activate picrust2v242 #( to use picrust2 v2.4.2 )

nohup picrust2_pipeline.py -s dada-chem.seqs.nochim2.fasta -i dada-chem.rcounts.ASV.2.tab --min_samples 0 -o picrust2_out_custom_tree -p 10 --ref_dir custom_ref/ --marker_gene_table trait.files.new/16S_SB.txt --no_pathways --skip_minpath --edge_exponent 0 &

# Unzip file
gunzip marker_predicted_and_nsti.tsv.gz
# Renaming it for downstream analysis in R
mv marker_predicted_and_nsti.tsv pic.16S_predicted_custom_tree.txt

