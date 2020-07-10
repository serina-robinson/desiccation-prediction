# Load necessary packages
# install.packages("pacman") # run if pacman not already installed
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("Biostrings")
pacman::p_load("Biostrings", "tidyverse") # if this throws an error, install packages first

# Set working directory
# setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/desiccation-prediction/") 
# You will change this to your own working directory, starting with ("C:/Users/...) 

# Read in example FASTA files
fasta1 <- readAAStringSet("data/test10.faa") # 10 sequences
fasta2 <- readAAStringSet("data/test5.faa") # 5 sequences
fasta3 <- readAAStringSet("data/test10.faa") # 10 sequences

# Combine all the FASTA files into one FASTA file
fasta_all <- AAStringSet(c(fasta1, fasta2, fasta3))
length(fasta_all) # Length is now 25 seqs (test10 + test5 + test10)

# Check if there are any duplicates
table(duplicated(fasta_all)) # there are 10 duplicates and 15 unique sequences
fasta_dedup <- fasta_all[!duplicated(fasta_all)] # here is how to remove the duplicated sequences
length(fasta_dedup) # now we have only 15 seqs (because test10 seqs were duplicated, 10 were removed)

# Now we are going to use these sequences as a query to makeblastdb and blastp 

# Let's pretend we are going to BLAST against these 119 genomes in the file data/nocardia_assemblies.txt
# You will use your own genome accessions!

# Read in the genome accessions
noc <- read.table('data/nocardia_assemblies.txt', stringsAsFactors = F)

# This is the basic syntax for each of the commands you'll be running for each genome

# On MSI:

# ssh yourusername@login.msi.umn.edu

# After password prompt, you will need to connect to mesabi:

# ssh mesabi
# qsub -l nodes=1:ppn=4,mem=4gb,walltime=4:00:00 -q interactive -I

# Once connected:
# module load ncbi_blast+

# curl all your Protein fasta files for your genomes of interest 
# This is just an example
# curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/145/GCF_000009145.1_ASM914v1_protein.faa.gz

# Then you will need to unzip them
# gunzip GCF_000009145.1_ASM914v1_protein.faa.gz


# For each fasta file make a blast database 
# Example:
# makeblastdb -in GCF_001613165.1_ASM161316v1_protein.faa -dbtype prot -out nocardia_db/GCF_001613165.1_ASM161316v1_protein_db

# Then run a Blastp with your query sequences
# Example:
# blastp -query query/20_query_seqs.fasta -db nocardia_db/GCF_001613165.1_ASM161316v1_protein_db -out blast_output/GCF_001613165.out -num_threads 4 -outfmt 6 -max_target_seqs 10

# Here is an example of how to write commmands for multiple genomes in a loop
# You can also do this for the curl and gunzip commands in a similar way

db <- list()
ll <- list()

for(i in 1:nrow(noc)){
  nocnam <- gsub("\\.faa",".db",noc[i,])
  db[[i]] <- paste0("makeblastdb -in ", noc[i,], " -dbtype prot -out nocardia_db/", nocnam)
  ll[[i]] <- paste0("blastp -query query/20_query_seqs.fasta -db nocardia_db/", nocnam, " -out blast_output/", nocnam, ".out -num_threads 4 -outfmt 6 -max_target_seqs 10")
}

unldb <- unlist(db)
unl <- unlist(ll)

# Write commands to text
write.table(unldb, "output/makeblastdb_commands.txt", row.names = F, quote = F, col.names = F)
write.table(unl, "output/nocardia_blastp_commands.txt", row.names = F, quote = F, col.names = F)

# And then transfer the files to MSI using Filezilla (barf)
# Then remember you have to use chmod u+x to make them executable files!
# For example:
# chmod u+x output/makeblastdb_commands.txt
# ./makeblastdb_commands.txt # should now run as a script!
