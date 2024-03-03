# This is a python script - works in conjunction with the 3. Compare Expected and Observed Homozygotes

import sys
import os
import pandas as pd
from pandas import DataFrame
import numpy as np

#sys.argv[1]=path to working directory 
#sys.argv[2]=file containing haplotype information e.g.joined.chr26.haplotype_block_0867.txt
#sys.argv[3]=file to write out e.g. homohapinfo.chr26.haplotype_block_0867.txt

#adjust this variable according to population size
pop_allele_count = int(85788)

#set up functions required

def find_homozygous_haplotypes(row):
 if (row['Haplotype1'] == row['Haplotype2']):
  return 'Homozygous'
 else:
  return 'Heterozygous'


#set working directory
os.chdir(sys.argv[1])

#read in file
haplotypes = pd.read_table(sys.argv[2], delimiter=' ', dtype = str)

#run function
haplotypes = haplotypes.assign(Diplotype=haplotypes.apply(find_homozygous_haplotypes, axis=1))
haplotypes = haplotypes[haplotypes['Diplotype'] == 'Homozygous']

#create table with frequencies of haplotypes
frequencies = haplotypes['Haplotype1'].value_counts().rename_axis('Haplotype').reset_index(name='Number_of_homo_haplotypes')
frequencies = frequencies.assign(total_homo_haplotype_observations=frequencies.Number_of_homo_haplotypes.mul(2))
frequencies = frequencies.assign(Observed_homozygous_frequency=frequencies.total_homo_haplotype_observations.div(pop_allele_count))
remove_counts = frequencies.pop('total_homo_haplotype_observations')

#write out file
frequencies.to_csv(sys.argv[3], sep='\t', index=False)
