import os
import subprocess

############################################
"""
Pair R1 and R2, align paired reads to genome using Bowtie2 and filtering for MAPQ of 23, make bam file, sort, index
    1. Looks for R1 and R2 .fastq.gz files in separate directories and with only string difference of "_R1_" and "_R2_" 
    in filename
    2. Makes "aligned" output directory with sub-folder "paired_reads"
    3. Seqkit pairs reads, dropping any unpaired that were individually removed from separate read trimming steps
    4. Performs paired alignment with paired R1 and R2 files, writing sam files to output directory
    5. Use samtools to convert sam to bam files and sort, create index
"""
############################################

#Location of bowtie2 reference genome
grch38 = '/seqdb/bowtie2/hg38'


#Define directories
R2_input_directory = "plasmid_filt"
R1_input_directory = "trimmed_fastq_R1"
output_directory = "aligned"
output_paired_reads_directory = os.path.join(output_directory, "paired_reads")

os.makedirs(output_directory, exist_ok=True)
os.makedirs(output_paired_reads_directory, exist_ok=True)


#Get a list of read1 and read2 FASTQ files in the input directories
read2_files = [filename for filename in os.listdir(R2_input_directory) if filename.endswith('.fastq.gz') and "_R2_" in filename]
read1_files = [filename for filename in os.listdir(R1_input_directory) if filename.endswith('.fastq.gz') and "_R1_" in filename]


#Pair and align reads using Bowtie2
def align_reads(read1_filename, read2_filename):  
    read2_file_path = os.path.join(R2_input_directory, read2_filename)
    read1_file_path = os.path.join(R1_input_directory, read1_filename)      

    #Get the base name of the files
    base = os.path.splitext(read1_filename.replace('_R1_', '_').replace('.fastq.gz', ''))[0]
    
    #Pair reads using seqkit
    seqkit_command = ["seqkit", "pair", 
    	"-1", read1_file_path, 
    	"-2", read2_file_path, 
    	"-O", output_paired_reads_directory,
    	"-j", "16"]

    try:
        subprocess.run(seqkit_command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running seqkit pair: {e}")
        return

    paired_read2_file_path = os.path.join(output_paired_reads_directory, read2_filename)
    paired_read1_file_path = os.path.join(output_paired_reads_directory, read1_filename)

    #Check if paired files are created
    if not os.path.exists(paired_read2_file_path) or not os.path.exists(paired_read1_file_path):
        print(f"Error: Paired files for {base} were not created.")
        return

    #Map the reads to the hg38 genome using Bowtie2
    bowtie2_command = ["bowtie2", "-x", grch38, "-1", paired_read1_file_path, "-2", paired_read2_file_path,
                       "--very-sensitive-local", "--threads", "16", "--no-unal", "-S",
                       os.path.join(output_directory, f"{base}.sam")]
    
    with open(os.path.join(output_directory, f"{base}_bowtieQC.txt"), "a") as f:
        subprocess.run(bowtie2_command, stdout=subprocess.PIPE, stderr=f)

#Check if both read files associated with a read are present
for read2_filename in read2_files:
    read1_filename = read2_filename.replace("_R2_", "_R1_")
    if read1_filename in read1_files:
        align_reads(read1_filename, read2_filename)
    else:
        print(f"{read2_filename}: file is missing read1 file")

print("Alignment complete")


####################################### SAM TO BAM ########################################

#Convert SAM files to BAM files and sort them
sam_files = [f for f in os.listdir(output_directory) if f.endswith('.sam')]
for file in sam_files:
    base_name = file.replace('.sam', '')
    bam_file_sorted_name = os.path.join(output_directory, f"{base_name}_sorted_name.bam")
    bam_file_sorted_coord = os.path.join(output_directory, f"{base_name}_sorted_coord.bam")

    #Convert SAM to BAM and sort by read names
    with open(bam_file_sorted_name, "wb") as output_file:
        process1 = subprocess.Popen(["samtools", "view", "-@", "16", "-bh", "-q", "23", os.path.join(output_directory, file)], stdout=subprocess.PIPE)
        process2 = subprocess.Popen(["samtools", "sort", "-@", "16", "-n", "-o", "-", "-"], stdin=process1.stdout, stdout=output_file)
        process1.stdout.close()
        process2.communicate()

    if process1.returncode != 0 or process2.returncode != 0:
        print(f"Error during SAM to BAM conversion for {base_name}")
        continue

    #Sort by genomic coordinates
    with open(bam_file_sorted_coord, "wb") as output_file:
        process3 = subprocess.Popen(["samtools", "sort", "-@", "16", "-o", "-", bam_file_sorted_name], stdout=output_file)
        process3.communicate()

    if process3.returncode != 0:
        print(f"Error sorting by coordinate for {base_name}")
        continue

    #Index the BAM file sorted by genomic coordinates
    subprocess.run(["samtools", "index", "-@", "16", bam_file_sorted_coord])

print("BAM files created and indexed")