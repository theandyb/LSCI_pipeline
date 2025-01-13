"""
Compare distributions of nucleotides flanking singletons across chromosomes
"""
import pandas as pd
import statsmodels.api as sm
#from pyfaidx import Fasta
import sqlite3
from patsy import dmatrices
import statsmodels.formula.api as smf
import argparse
import time

parser = argparse.ArgumentParser(description="Compare flanking nucleotide distributions")
parser.add_argument("-p", "--population", help="1kGP Superpopulation Code", required=True)
parser.add_argument("-s", "--subtype", help="Mutation Subtype", required=True)
parser.add_argument("-o", "--output", help="Output base directory", required=True)
parser.add_argument("-x", "--suffix", help="Control file suffix", default="")
parser.add_argument("-c", "--chrom", help="First chromosome", default="")
parser.add_argument("-d", "--drom", help="Second chromosome", default="")
args = parser.parse_args()

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

def get_seq_db(chromosome, db_file="data/reference/ref_db.db"):
  db_file_ro = "file:{}?mode=ro".format(db_file)
  conn = sqlite3.connect(db_file_ro, uri=True)
  c = conn.cursor()
  sequence_id = f"chr{chromosome}"
  c.execute('SELECT sequence FROM sequences WHERE id=?', (sequence_id,))
  result = c.fetchone()
  conn.close()
  return(result[0])

def get_count_table_singletons(chromosome, subtype, offset, pop, base_dir):
  """
  Get the count table for singletons at a given relative position for a subtype.
  Note: this version has been made so as to be called in parallel via ray, but
  that it sucks in that each worker has to load the reference genome. Need to
  check if there's a python library that has the reference genome as on object
  that can maybe be passed around?
  """
  singleton_pos = s_pos(subtype, chromosome, pop, base_dir)
  rev_singleton_pos = s_pos(subtype + "_rev", chromosome, pop, base_dir)
  #fasta_obj = Fasta(ref_file)
  #seq = fasta_obj["{}{}".format("chr", chromosome)]
  #seqstr = seq[0:len(seq)].seq
  seq = get_seq_db(chromosome)
  results = {"A":0, "C":0, "G":0, "T":0}
  for index, value in singleton_pos.items():
    ix = value - 1 + offset
    nuc = seq[ix]
    if nuc in results.keys():
      results[nuc] += 1
  for index, value in rev_singleton_pos.items():
    ix = value - 1 + (offset * -1)
    nuc = complement(seq[ix])
    if nuc in results.keys():
      results[nuc] += 1
  return(results)

def fit_model(chrom1, chrom2, subtype, offset, pop, base_dir, suffix = ""):
  s_tab = get_count_table_singletons(chrom1, subtype, offset, pop, base_dir)
  s_tab = pd.DataFrame(s_tab.items(), columns=['nuc', 'singletons'])
  c_tab = get_count_table_singletons(chrom2, subtype, offset, pop, base_dir)
  c_tab = pd.DataFrame(c_tab.items(), columns=['nuc', 'controls'])
  df = pd.DataFrame.merge(s_tab, c_tab, on='nuc')
  df = pd.melt(df, id_vars = 'nuc', var_name = 'status', value_name = "n")
  mod = smf.glm("n ~ status + nuc", df, family = sm.families.Poisson()).fit()
  df['fitted'] = mod.fittedvalues
  df['res'] = mod.resid_deviance
  return(df)

population = args.population
subtype = args.subtype
base_dir = args.output
suffix = args.suffix
chromosome1 = str(args.chrom)
chromosome2 = str(args.drom)

print("Running with options: ")
print("Population: " + population)
print("Subtype: " + subtype)
print("Base directory: " + base_dir)
print("Chromosome 1: " + chromosome1)
print("Chromosome 2: " + chromosome2)


res_out_dir = "{}/single_pos/chrom_comp/".format(base_dir)
print("Running models for subtype: {} and population: {}".format(subtype, population))
for offset in range(1, 21):
  print(offset)
  df = fit_model(chromosome1, chromosome2, subtype, offset * -1, population, base_dir, suffix = suffix)
  #fit_model_all(subtype, offset, pop, ref_file, base_dir, suffix = "")
  file_name = "{}{}_{}_{}_{}_{}.csv".format(res_out_dir, population, subtype, chromosome1, chromosome2, str(-1*offset))
  df.to_csv(file_name, index = False)
  if subtype.startswith("cpg") and offset == 1:
    continue
  df = fit_model(chromosome1, chromosome2, subtype, offset, population, base_dir, suffix = suffix)
  file_name = "{}{}_{}_{}_{}_{}.csv".format(res_out_dir, population, subtype, chromosome1, chromosome2, str(offset))
  df.to_csv(file_name, index = False)
