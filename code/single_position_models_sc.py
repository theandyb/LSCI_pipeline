"""
Single Chromosome Results
Code for getting the deviance statistics at all positions within +/- 1000bp window
"""
import pandas as pd
import statsmodels.api as sm
from pyfaidx import Fasta
from patsy import dmatrices
import statsmodels.formula.api as smf
import argparse
import time

parser = argparse.ArgumentParser(description="Fit single position models")
parser.add_argument("-p", "--population", help="1kGP Superpopulation Code", required=True)
parser.add_argument("-s", "--subtype", help="Mutation Subtype", required=True)
parser.add_argument("-f", "--fasta", help="Reference Genome Fasta", required=True)
parser.add_argument("-o", "--output", help="Output base directory", required=True)
parser.add_argument("-x", "--suffix", help="Control file suffix", default="")
parser.add_argument("-c", "--chrom", help="Chromosome", default="")
args = parser.parse_args()

def c_pos(subtype, chromosome, pop, base_dir, suffix):
  """Get the positions for controls for a given subtype"""
  input_dir = "{}/controls/{}/pos_files/".format(base_dir, pop)
  f_name = input_dir + subtype + "_" + str(chromosome) + ".txt" + suffix
  pos_list = pd.read_csv(f_name, header=None, names = ['pos'], usecols=['pos']).squeeze("columns")
  return pos_list

def s_pos(subtype, chromosome, pop, base_dir):
  """Get the positions for singletons for a given subtype"""
  input_dir = "{}/singletons/{}/pos_files/".format(base_dir, pop)
  f_name = input_dir + subtype + "_" + str(chromosome) + ".txt"
  pos_list = pd.read_csv(f_name, header=None, names = ['pos'], usecols=['pos']).squeeze("columns")
  return pos_list

def complement(nucleotide):
  complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
  if not nucleotide in complements.keys():
    return nucleotide
  return complements[nucleotide]

def get_count_table_control(chromosome, subtype, offset, pop, ref_file, base_dir, suffix):
  """
  Get the count table for controls at a given relative position for a subtype.
  Note: this version has been made so as to be called in parallel via ray, but
  that it sucks in that each worker has to load the reference genome. Need to
  check if there's a python library that has the reference genome as on object
  that can maybe be passed around?
  """
  fasta_obj = Fasta(ref_file)
  control_pos = c_pos(subtype, chromosome, pop, base_dir, suffix = suffix)
  rev_control_pos = c_pos(subtype + "_rev", chromosome, pop, base_dir, suffix = suffix)
  seq = fasta_obj["{}{}".format("chr", chromosome)]
  #seqstr = seq[0:len(seq)].seq
  results = {"A":0, "C":0, "G":0, "T":0}
  for index, value in control_pos.items():
    ix = value - 1 + offset
    nuc = seq[ix].seq
    if nuc in results.keys():
      results[nuc] += 1
  for index, value in rev_control_pos.items():
    ix = value - 1 + (offset * -1)
    nuc = complement(seq[ix].seq)
    if nuc in results.keys():
      results[nuc] += 1
  return(results)

def get_count_table_singletons(chromosome, subtype, offset, pop, ref_file, base_dir):
  """
  Get the count table for singletons at a given relative position for a subtype.
  Note: this version has been made so as to be called in parallel via ray, but
  that it sucks in that each worker has to load the reference genome. Need to
  check if there's a python library that has the reference genome as on object
  that can maybe be passed around?
  """
  singleton_pos = s_pos(subtype, chromosome, pop, base_dir)
  rev_singleton_pos = s_pos(subtype + "_rev", chromosome, pop, base_dir)
  fasta_obj = Fasta(ref_file)
  seq = fasta_obj["{}{}".format("chr", chromosome)]
  #seqstr = seq[0:len(seq)].seq
  results = {"A":0, "C":0, "G":0, "T":0}
  for index, value in singleton_pos.items():
    ix = value - 1 + offset
    nuc = seq[ix].seq
    if nuc in results.keys():
      results[nuc] += 1
  for index, value in rev_singleton_pos.items():
    ix = value - 1 + (offset * -1)
    nuc = complement(seq[ix].seq)
    if nuc in results.keys():
      results[nuc] += 1
  return(results)


def fit_model(chromosome, subtype, offset, pop, ref_file, base_dir, suffix = ""):
  s_tab = get_count_table_singletons(chromosome, subtype, offset, pop, ref_file, base_dir)
  s_tab = pd.DataFrame(s_tab.items(), columns=['nuc', 'singletons'])
  c_tab = get_count_table_control(chromosome, subtype, offset, pop, ref_file, base_dir, suffix)
  c_tab = pd.DataFrame(c_tab.items(), columns=['nuc', 'controls'])
  df = pd.DataFrame.merge(s_tab, c_tab, on='nuc')
  df = pd.melt(df, id_vars = 'nuc', var_name = 'status', value_name = "n")
  mod = smf.glm("n ~ status + nuc", df, family = sm.families.Poisson()).fit()
  df['fitted'] = mod.fittedvalues
  df['res'] = mod.resid_deviance
  return(df)

ref_genome = args.fasta
population = args.population
subtype = args.subtype
base_dir = args.output
suffix = args.suffix
chromosome = str(args.chrom)

print("Running with options: ")
print("Population: " + population)
print("Subtype: " + subtype)
print("Reference File: {}".format(ref_genome))
print("Base directory: " + base_dir)
print("Chromosome: " + chromosome)

time.sleep(3)

res_out_dir = "{}/single_pos/chr{}/resid/{}/".format(base_dir, chromosome, population)
print("Running models for subtype: {} and population: {}".format(subtype, population))
for offset in range(1, 501):
  print(offset)
  df = fit_model(chromosome, subtype, offset * -1, population, ref_genome, base_dir, suffix = suffix)
  #fit_model_all(subtype, offset, pop, ref_file, base_dir, suffix = "")
  file_name = res_out_dir + subtype + "_rp_" + str(-1*offset) + ".csv" + suffix
  df.to_csv(file_name, index = False)
  if subtype.startswith("cpg") and offset == 1:
    continue
  df = fit_model(chromosome, subtype, offset, population, ref_genome, base_dir, suffix = suffix)
  file_name = res_out_dir + subtype + "_rp_" + str(offset) + ".csv" + suffix
  df.to_csv(file_name, index = False)

