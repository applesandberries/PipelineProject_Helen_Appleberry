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
