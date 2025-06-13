#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 17:54:08 2023

@author: smelab
"""

import glob, re, os
import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt
from Bio import SeqIO, Entrez
from Bio.Seq import Seq

def gene_informn(entrez_id):
    '''
    This code gets information for a given entrez id from NCBI database. Stable
    internet connection is needed for this code

    Parameters
    ----------
    entrez_id : String
        The entrez id for the gene.

    Returns
    -------
    gene_info : Dictionary
        Return Start codon pos, stop codon pos, second in-frame stop codon pos
        and transcript length for the input entrez id. Uncomment the lines for
        other details if needed

    '''
    Entrez.email = 'saubhiksom@iisc.ac.in'
    pattern_1 = "\d+"
    with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=entrez_id) as handle:
        seq_record = SeqIO.read(handle, "gb") # using "gb" as an alias for "genbank"
    gene_info = {}               
    for features in seq_record.features:
        if features.type == "CDS":
            transcript_seq = str(seq_record.seq)
            location = str(features.location)
            CDS_location_ATG = int(re.findall(pattern_1, location)[0])
            CDS_location_STOP = int(re.findall(pattern_1, location)[1])
            UTR_seq = str(seq_record.seq[CDS_location_STOP:])
            ISR_protein = str(Seq.translate(seq_record.seq[CDS_location_STOP:], to_stop=True))
            gene_info[entrez_id] = {'start': CDS_location_ATG,
                                   "stop1": CDS_location_STOP,
                                   "stop2": CDS_location_STOP+len(ISR_protein)*3+3,
                                   "length_gene": len(transcript_seq),
                                   "UTR_length": len(UTR_seq)}
            break
    return gene_info

def plotting_sorted_BAM(gene_name, entrez_id, list_of_files=None, g4_coords=None):
    '''
    This function counts the reads from a sorted BAM file and returns a plot.

    Parameters
    ----------
    gene_name : Name of the gene. A folder will be created in the working 
                directory with same name
    entrez_id : The NM_xxxxx id of the gene
    
    list_of_files: Give list of SRRs if you want. Those files must be in sorted
                    BAM file folder
    
    g4_coords: Dictionary with G4 coordinates for each gene

    Returns:
        Plots the numbers of reads per million for the input transcript
        returns a pandas dataframe with alignment info in CDS and UTR
    -------
    '''
    
    gene_info = gene_informn(entrez_id)
    g4_start, g4_end = g4_coords[entrez_id]
    start = gene_info[entrez_id]['start']
    stop1 = gene_info[entrez_id]['stop1']
    stop2 = gene_info[entrez_id]['stop2']
    length_gene = gene_info[entrez_id]["length_gene"]

    # Create output location for plots
    location = f"{gene_name}_plots/"
    if os.path.isdir(location) == False:
        os.makedirs(location)
    
    # Look for BAM files in sorted_bam_files folder
    bam_files = glob.glob("sorted_bam_files/*.sorted.bam")
    
    if list_of_files != None:
        if type(list_of_files) != list:
            print("List input needed")
            return
        else:
            bam_files = [i for i in bam_files if (re.findall('SRR\d*', i)[0]) in list_of_files]
    
    if len(bam_files) == 0:
        print(f"No bam files found in sorted_bam_files/")
        return
    
    total_reads_all_srrs = 0
    a_total = np.zeros(length_gene)
    all_tables = None
    
    for filename in bam_files:
        srr = re.findall('SRR\d*', filename)[0]
        a = np.zeros(length_gene)
        sorted_bamfile = pysam.AlignmentFile(filename, "rb")
        total_reads = 0
        
        for read in sorted_bamfile.fetch(until_eof=True):
            total_reads += 1
            total_reads_all_srrs += 1
            try:
                if read.reference_name == entrez_id:
                    position = read.pos
                    length = read.qlen
                    a[position - 1 : position - 1 + length] += 1
                    a_total[position - 1 : position - 1 + length] += 1
            except ValueError:
                continue

        sorted_bamfile.close()        
        
        a_rpm = (a / total_reads) * 1e6 #calculating RPM
        fig, ax = plt.subplots(1, 1, figsize=[16, 8], dpi=600)
        ax.plot(a_rpm)
        ax.set_xticks(ticks=[start, stop1, stop2])
        ax.set_xticklabels(labels=['S', '*', '*'])
        ax.set_ylabel('Reads per million', fontsize=18)
        ax.set_xlabel("Codon position", fontsize=18)
        fig.suptitle(gene_name + '_' + srr, fontsize=20)
        
        # Save plot in current directory
        plot_filename = f"{location}{gene_name}_{srr}.png"
        fig.savefig(plot_filename, dpi=300)
        plt.show()
        
        # Calculate G4 region statistics
        up_start = max(0, g4_start - 30)
        up_end = g4_start - 1
        down_start = g4_end + 1
        down_end = min(length_gene - 1, g4_end + 30)
        rd_g4 = np.sum(a[g4_start:g4_end+1]) / (g4_end - g4_start + 1)
        rd_up = np.sum(a[up_start:up_end+1]) / (up_end - up_start + 1)
        rd_down = np.sum(a[down_start:down_end+1]) / (down_end - down_start + 1)
        
        table_data = pd.DataFrame([[entrez_id, gene_name, srr, length_gene, rd_up, rd_g4, rd_down]], 
                                 columns=["Gene ID", "Gene Name", "SRR No", "Length", "RPBM Upstream", "RPBM G4", "RPBM Downstream"])
        
        if all_tables is None:
            all_tables = table_data
        else:
            all_tables = pd.concat([all_tables, table_data], ignore_index=True)
    
    return all_tables
