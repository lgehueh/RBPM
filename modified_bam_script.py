#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 11:04:33 2023

@author: smelab
"""

import os
import pysam
import subprocess
import pandas as pd
from plotting_reads_from_sorted_BAM_file import plotting_sorted_BAM

# Set your working directory - ALL outputs will go here
WORKING_DIR = "/media/smelab/New_Volume/kritika"
os.chdir(WORKING_DIR)

# Define the SRA accession numbers and transcript filename
gene_name = 'sorted_bam_files'

sra_accessions = ['SRR970587', 'SRR970538','SRR970565','SRR970490','SRR970588','SRR9113069','SRR9113063','SRR403885', 'SRR1248253', 'SRR1598971', 'SRR1573934', 'SRR2064017','SRR2064020','SRR4450327']

# Transcript file is in the same directory (from git repo)
transcript_filename = "transcript11.fasta"

g4_coords = {
    "NM_004565.3":(1047,1076),
    "NM_000418.4":(2057,2086),
    "NM_014727.3":(288,309),
    "NM_014727.3":(3157,3185),
    "NM_001164586.2":(7837,7866),
    "NM_001033910.3":(1407,1424),
    "NM_001077268.2":(614,636),
    "NM_020832.3":(1467,1496),
    "NM_052857.4":(522,551),
    "NM_015355.4":(351,380),
    "NM_024619.4":(419,448),
    "NM_002024.6":(1810,1838)
}

genes_to_process = [
    ("PEX14","NM_004565.3"),
    ("IL4R","NM_000418.4"),
    ("MLL4","NM_014727.3"),
    ("IGFN1","NM_001164586.2"),
    ("TRAF5","NM_001033910.3"),
    ("ZFYVE19","NM_001077268.2"),
    ("ZNF687","NM_020832.3"),
    ("ZNF830","NM_052857.4"),
    ("SUZ12","NM_015355.4"),
    ("FN3KRP","NM_024619.4"),
    ("FMR1","NM_002024.6")
]

# Build the bowtie2 index for the transcript (in working directory)
transcript_index = "transcript"
bowtie2_build_command = ["/home/smelab/bowtie2/bowtie2-build", transcript_filename, transcript_index]
subprocess.check_call(bowtie2_build_command)
print("Bowtie2 index building done")

# Download SRA and Align the reads to the transcript
aligned_reads = []
total_reads = 0

for sra_accession in sra_accessions:
    # All files will be in the working directory
    fastq_filename = f"{sra_accession}.fastq"
    
    if os.path.isfile(fastq_filename) == False:
        print(f'Starting download for {sra_accession}')
        fastq_dump_command = ["/home/smelab/sratoolkit.3.0.5-ubuntu64/bin/fastq-dump", sra_accession, '--outdir', '.']
        subprocess.check_call(fastq_dump_command)
        print('Download fastq done')

    sam_filename = f"{sra_accession}_aligned_reads.sam"
    bowtie2_command = ["/home/smelab/bowtie2/bowtie2", "-x", 
                       transcript_index, "-U", 
                       fastq_filename, "-S", 
                       sam_filename, "-p", "1", "--very-sensitive-local"]
    try:
        subprocess.run(bowtie2_command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running Bowtie 2: {e}")
        continue
    
    # Convert the SAM file to a BAM file
    bam_filename = f"{sra_accession}_aligned_reads.bam"
    samfile = pysam.AlignmentFile(sam_filename, "r")
    bamfile = pysam.AlignmentFile(bam_filename, "wb", template=samfile)
    for read in samfile:
        bamfile.write(read)
    samfile.close()
    bamfile.close()
    
    # Sort the BAM file in sorted_bam_files subdirectory
    sorted_bam_dir = gene_name
    sorted_bam_filename = f"{sorted_bam_dir}/{sra_accession}_aligned_reads.sorted.bam"
    os.makedirs(sorted_bam_dir, exist_ok=True)

    pysam.sort("-o", sorted_bam_filename, bam_filename)
    os.remove(bam_filename)

    # Index the sorted BAM file
    pysam.index(sorted_bam_filename)

    # Delete intermediate files
    os.remove(sam_filename)
      
    # Get alignment information
    aligned_reads.append((sra_accession, sorted_bam_filename))

    # Print the number of mapped reads for each SRA
    sorted_bamfile = pysam.AlignmentFile(sorted_bam_filename, "rb")
    num_mapped_reads = 0
    for read in sorted_bamfile:
        if not read.is_unmapped:
            num_mapped_reads += 1
    sorted_bamfile.close()
    print(f"Number of mapped reads for {sra_accession}: {num_mapped_reads}")
    total_reads += num_mapped_reads

# Print the total number of mapped reads across all SRAs
print(f"Total number of mapped reads: {total_reads}")

# Process genes and generate Excel output
all_gene_data = []

for gene_name, entrez_id in genes_to_process:
    print(f"\nProcessing {gene_name} ({entrez_id})...")
    if entrez_id not in g4_coords:
        print(f"G4 coordinates not defined for {entrez_id}. Skipping.")
        continue
    df = plotting_sorted_BAM(
        gene_name,
        entrez_id,
        list_of_files=sra_accessions,
        g4_coords=g4_coords  # Pass G4 coordinates
    )
    all_gene_data.append(df)

# Merge and save all results in working directory
if all_gene_data:
    final_df = pd.concat(all_gene_data, ignore_index=True)
    output_excel = "G4_RPBM_All_Genes.xlsx"
    final_df.to_excel(output_excel, index=False)
    print(f"Saved: {output_excel}")
else:
    print("No gene data was processed.")

# Clean up transcript index files from working directory
index_files = [filename for filename in os.listdir('.') if filename.startswith("transcript.")]
for index_file in index_files:
    os.remove(index_file)

print(f"Cleaned up {len(index_files)} index files")