import subprocess
import os
import gzip
from Bio import Entrez
from Bio import SeqIO
import glob
# Get the current working directory
basic_path = os.getcwd()
Entrez.email = 'helen.m.appleberry@gmail.com'
# Define the file path for the fastq files and the other files from the ncbi fasta
file_path=basic_path+'/fastq_files/'
more_files=basic_path+'/more_files'
#doing a little blast for the nucleotide provided
handle = Entrez.efetch(db="nucleotide", id= 'NC_006273.2', rettype='fasta')
records = list(SeqIO.parse(handle, 'fasta'))
SeqIO.write(records, more_files, "fasta")
log=open(basic_path+'/PipelineProject_Helen_Appleberry/PipelineProject.log', 'w+')
#os.makedirs(basic_path+'/index')
index = basic_path+'/index'
os.system('bowtie2-build ' + more_files + " "+ index)
index = basic_path+'/index'
#lets grab all the files!
x=glob.glob (file_path+'*'+'.fastq')
list_of_fastq_prefixes=set()
#making a list for all the prefixes
for i in x:
    p=i[len(file_path):i.find('_1')]
    if '_2' not in i and '.fastq' in i:
        list_of_fastq_prefixes.add(p)
print(list_of_fastq_prefixes)
#doing a bowtie command for each file 1 and 2 so we get a mapped!
for i in list_of_fastq_prefixes:
    bowtie2_command = f'bowtie2 -x {index} -1 {os.path.join(file_path, f"{i}_1.fastq")} -2 {os.path.join(file_path, f"{i}_2.fastq")} -S {i}_map.sam --al-conc {i}_mapped_%.fq'
    print(bowtie2_command)
    #os.system(bowtie2_command)
#okay so now we need to open the fastq files before mapping and count the number of reads so that we can write that amount to our log file
for i in list_of_fastq_prefixes:
    #Count read pairs before mapping
    with open(os.path.join(file_path, f"{i}_1.fastq"), "r") as fastq_file:
        records_before = list(SeqIO.parse(fastq_file, "fastq"))
        total_read_pairs_before = (len(records_before))
# Count read pairs after mapping
    with open(os.path.join(basic_path, f"{i}_mapped_1.fq"), "r") as file:
        file=file.readlines()
        total_read_pairs_after =(len(file) // 4)
#write this to the log file
    log.write(str(f'Donor {i} had {total_read_pairs_before} read pairs before Bowtie2 filtering and {total_read_pairs_after} read pairs after mapping'))
#okay so at this point I realized my prefixes needed to be a list in order to iterate through in the next step
list_of_fastq_prefixes = list(list_of_fastq_prefixes)
#making a path for where I am writing these assembly files
path = basic_path + '/Assemblies/'

#my amazing spades command where in I iterate through my list of prefixes and make lovely little mapped files

spades_command = f'rnaspades.py -t 2  --pe-1 1 '+basic_path+'/'+list_of_fastq_prefixes[0]+'_mapped_1.fq  --pe-2 1 '+basic_path+'/'+list_of_fastq_prefixes[0]+'_mapped_2.fq  --pe-1 2 '+basic_path+'/'+list_of_fastq_prefixes[1]+'_mapped_1.fq --pe-2 2 '+basic_path+'/'+list_of_fastq_prefixes[1]+'_mapped_2.fq --pe-1 3 '+basic_path+'/'+list_of_fastq_prefixes[2]+'_mapped_1.fq --pe-2 3 '+basic_path+'/'+list_of_fastq_prefixes[2]+'_mapped_2.fq --pe-1 4 '+basic_path+'/'+list_of_fastq_prefixes[3]+'_mapped_1.fq --pe-2 4 '+basic_path+'/'+list_of_fastq_prefixes[3]+'_mapped_2.fq -o assembly1'

log.write(spades_command)
#os.system(spades_command)
#makinga path to the transcripts
assembly_file= basic_path+ '/assembly1/transcripts.fasta'
#initializing my count
count=0
#opening the assembly file and counting the number that are larger than 1000
with open(assembly_file, "r") as assembly_handle:
    for record in SeqIO.parse(assembly_handle, "fasta"):
        if len(record.seq) > 1000:
            count += 1
log.write(str(f'There are {count} contigs > 1000 bp in the assembly.'))
#initializing the total length veriable
total_length = 0
#counting those numbers of the files which are over 1000!
#
with open(assembly_file, "r") as assembly_handle:
    for record in SeqIO.parse(assembly_handle, "fasta"):
        if len(record.seq) > 1000:
            total_length += len(record.seq)

log.write(str(f'There are {total_length} bp in the assembly.'))
#this makes the folder for blast data
#os.makedirs(basic_path+'/blast')
#initializing our variables
max_len = 0
max_contig = ''

#this finds the longest contig in transcripts.fasta
with open(assembly_file, "r") as assembly_handle:
    for record in SeqIO.parse(assembly_handle, "fasta"):
        if len(record.seq) > max_len:
            max_len = len(record.seq)
            max_contig = str(record.seq)
#writes out longest contig into separate fasta file
output_fasta_path = os.path.join(basic_path, 'blast', 'longest_contig.fasta')
with open(output_fasta_path ,'w') as longest_contig:
    longest_contig.write(max_contig)
#making a path to the betaherpesvirinae
beta_db_path = os.path.join(basic_path, 'blast', 'beta_db.fasta')
#opening this file and writing the results of the blast to it.
with open(beta_db_path, 'w') as beta_db:
    handle = Entrez.esearch(db='nucleotide', term='Betaherpesvirinae[Organism] OR Betaherpesvirinae[All Fields]', retmax=1500)
    record = Entrez.read(handle)
    ids = record["IdList"]
    handle = Entrez.efetch(db='nucleotide', id=ids, rettype='fasta')

    records = list(SeqIO.parse(handle, format='fasta'))
#i needed it to all be in one fasta file so a looped it all through
    for record in records:
        beta_db.write(">" + str(record.description) + "\n")
        beta_db.write(str(record.seq) + "\n")
beta_db.close()
#i wrote the blast command with it called db_name so instead of correcting it i thought it would be better to just rename it?? Idk I did this part when i had a fever.
db_name=beta_db_path
#make the path for the mkblast and then make the blast
out_mk_blast= os.path.join(basic_path, 'blast', 'mkblast')
makeblast_command = f'makeblastdb -in '+db_name+' -out '+out_mk_blast+' -title '+out_mk_blast+' -dbtype nucl'
os.system(makeblast_command)
#making the blast output csv!!!
blast_output_csv = (basic_path+'/blast/blast_output.csv')
#again why am I moving the directory?? Dont question fever induced helen
os.chdir(os.path.join(basic_path,'blast'))
#DOing the last big command! a lot of googling went into realizing that you gave us what we need in the literal problem
blast_command = 'blastn -query '+output_fasta_path+' -db '+out_mk_blast+" -num_threads 2 -max_hsps 1 -max_target_seqs 10 -out "+blast_output_csv+" -outfmt='6 sacc pident length qstart qend sstart send bitscore evalue stitle'"
os.system(blast_command)
log.write(str("sacc" + "\t" + "pident" + "\t" + "length" + "\t" +  "qstart" + "\t" + "qend" + "\t" + "sstart" + "\t" + "send" + "\t" + "bitscore" + "\t" + "evalue" + "\t" + "stitle" + "\n"))
with open(blast_output_csv, 'r') as blast_file:
    content = blast_file.read()
    log.write(content)
log.close()
print(f"Blast results written to {blast_output_csv}")

