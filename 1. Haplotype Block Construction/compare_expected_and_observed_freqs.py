# This is a python script - works in conjunction with 3. Compare Expected and Observed Homozygotes

import os
import pandas as pd
from pandas import DataFrame
import numpy as np
import scipy
from scipy import stats
from scipy.stats import binom
import sys
from mpmath import *

#Inputs required:
#sys.argv[1]=work directory e.g. /data/seq/swjiv0/haplotype_analysis/Haplotype_Analysis_March2022/haplotype_blocks/chr26
#sys.argv[2]=file containing information on expected haplotype frequencies e.g. chr26.haplotype_frequencies.txt
#sys.argv[3]=file containing information on observed haplotype frequencies e.g. chr26.homozygous_haplotype_info.txt
#sys.argv[4]=output file e.g. chr26_candidate_hom_depleted_haplotypes.txt

#adjust this variable according to population size
pop_size = int(42894)

#set up functions required
def missing_homozygotes(row):
 if(float(row['Expected_number_of_homo_haplotypes']) >= float(10)) & (float(row['Number_of_homo_haplotypes']) == float(0)):
  return 'Missing_homozygotes'
 else:
  return 'Deficit_homozygotes'

def remove_common_haplotypes(row):
 if (float(row['Haplotype_frequency']) > float(0.15)):
  return 'Common'
 else:
  return 'Rare'

def compare_expected_vs_observed_homos(row):
 if (int(row['Expected_number_of_homo_haplotypes']) >= int(3)):
  if (int(row['Number_of_homo_haplotypes']) < int(row['Expected_number_of_homo_haplotypes'])):
   return 'Less homozygotes than expected'
  else:
   return 'More homozygotes than expected'
 else:
  return 'Not enough observations'

def homozygous_depletion(df):
 probability = []
 for index, row in df.iterrows():
  x = binom.pmf(int(row['Number_of_homo_haplotypes'] ) , pop_size, float(row['Expected_homozygous_frequency'] ))
  probability.append(x)
 df['Binomial_probability'] = probability

#NOTE: in the function called "homozygous_depletion" the number of expected observations = population size
#In this example population size = 42,894

#set working directory
os.chdir(sys.argv[1])

#read in file
hap_info = pd.read_table(sys.argv[2], delimiter='\t', dtype = {'HapBlockID' : str, 'Haplotype' : str, 'Count' : int, 'Haplotype_frequency' : float, 'Expected_homozygous_frequency' : float})
frequencies = pd.read_table(sys.argv[3], delimiter='\t', dtype={'Haplotype' : str, 'Number_of_homo_haplotypes' : int, 'Observed_homozygous_frequency' : float})

#Calculate the number of animals we would expect to be homozygous for each hapltoype
#NOTE the multiplier = population count. This may have to be changed based on your pop size
hap_info = hap_info.assign(Expected_number_of_homo_haplotypes=round(hap_info.Expected_homozygous_frequency.mul(pop_size), 0))

#merge haplotype info - expected and observed
homo_hap_freqs = pd.merge(hap_info, frequencies, on='Haplotype', how='outer')
homo_hap_freqs.fillna('0', inplace=True)

#identify which haplotypes are missing and which are depleted in their homozygous state
homo_hap_freqs = homo_hap_freqs.assign(Type=homo_hap_freqs.apply(missing_homozygotes, axis=1))

#remove haplotypes common to the population
homo_hap_freqs = homo_hap_freqs.assign(Frequency_filter=homo_hap_freqs.apply(remove_common_haplotypes, axis=1))
homo_hap_freqs = homo_hap_freqs[homo_hap_freqs['Frequency_filter'] == 'Rare']
remove_col=homo_hap_freqs.pop('Frequency_filter')

#only keep haplotypes that where there are less than expected homozygotes
homo_hap_freqs = homo_hap_freqs.assign(Compare_num_homos=homo_hap_freqs.apply(compare_expected_vs_observed_homos, axis=1))
homo_hap_freqs = homo_hap_freqs[homo_hap_freqs['Compare_num_homos'] == 'Less homozygotes than expected']
remove_col=homo_hap_freqs.pop('Compare_num_homos')

#use binomial probability function to see if the difference between expected and observed are statstically significant
homozygous_depletion(homo_hap_freqs)

#write out results
homo_hap_freqs.to_csv(sys.argv[4], sep='\t', index=False)
