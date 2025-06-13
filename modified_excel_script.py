#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 17:04:20 2023

@author: smelab
"""
import os
import pandas as pd
from plotting_reads_from_sorted_BAM_file import plotting_sorted_BAM

# Set your working directory - ALL outputs will go here
WORKING_DIR = "/media/smelab/New_Volume/kritika"
os.chdir(WORKING_DIR)

# Define G4 coordinates (if needed for your specific analysis)
g4_coords = {
    'NM_012199.5': (1000, 1030)  # Add your AGO1 G4 coordinates here
}

ago = ['SRR1605309','SRR2052945','SRR5013257','SRR2096965','SRR5345623','SRR4293695']
data = []

for sra in ago:
    result = plotting_sorted_BAM('AGO1', 'NM_012199.5', [sra], g4_coords=g4_coords)
    if result is not None:
        data.append(result)
    
if data:
    excel_data = pd.concat(data, axis=0, ignore_index=True)
    output_file = 'ago1_excel_data.xlsx'  # Will be saved in current directory
    excel_data.to_excel(output_file, index=False)
    print(f"Excel file saved: {output_file}")
else:
    print("No data to save")