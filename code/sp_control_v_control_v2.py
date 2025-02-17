"""
Code for getting the deviance statistics at all positions within +/- 1000bp window
"""
import pandas as pd
import statsmodels.api as sm
#from pyfaidx import Fasta
import sqlite3
from patsy import dmatrices
import statsmodels.formula.api as smf
import ray
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-p", "--population", help="1kGP Super Population", required=True, choices=["ALL","AFR","AMR","EAS","EUR","SAS"])
parser.add_argument("-t", "--subtype", help="Mutation subtype", choices=["AT_CG", "AT_GC", "AT_TA", "GC_AT", "GC_TA", "GC_CG", "cpg_GC_AT", "cpg_GC_TA", "cpg_GC_CG"], required=True)
parser.add_argument("-o", "--output", help="Output base directory", required=True)
args = parser.parse_args()

# Set these parameters before running
pop = args.population
subtype = args.subtype
base_dir = args.output

def c_pos(subtype, chromosome, pop, base_dir, start=0, skip=5):
  """Get the positions for controls for a given subtype"""
  input_dir = "{}/controls/{}/pos_files/".format(base_dir, pop)
  f_name = input_dir + subtype + "_" + str(chromosome) + ".txt"
  pos_list = pd.read_csv(f_name, header=None, names = ['pos'], usecols=['pos']).squeeze("columns")
  return pos_list[start::skip].astype(int)

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

@ray.remote
def get_count_table_control(chromosome, subtype, offset, pop, base_dir, start=0, skip =5):
  """
  Get the count table for controls at a given relative position for a subtype.
  Note: this version has been made so as to be called in parallel via ray, but
  that it sucks in that each worker has to load the reference genome. Need to
  check if there's a python library that has the reference genome as on object
  that can maybe be passed around?
  """
  #fasta_obj = Fasta(ref_file)
  control_pos = c_pos(subtype, chromosome, pop, base_dir, start = start, skip = skip)
  rev_control_pos = c_pos(subtype + "_rev", chromosome, pop, base_dir, start = start, skip = skip)
  #seq = fasta_obj["{}{}".format("chr", chromosome)]
  #seqstr = seq[0:len(seq)].seq
  seq = get_seq_db(chromosome)
  results = {"A":0, "C":0, "G":0, "T":0}
  for index, value in control_pos.items():
    ix = value - 1 + offset
    nuc = seq[ix]
    if nuc in results.keys():
      results[nuc] += 1
  for index, value in rev_control_pos.items():
    ix = value - 1 + (offset * -1)
    nuc = complement(seq[ix])
    if nuc in results.keys():
      results[nuc] += 1
  return(results)

def get_count_all(subtype, offset, pop, base_dir,  status = "singleton"):
  """
  Get the nearest and furthest control counts for a given relative position
  across all 22 autosomes. This version is parallelized via ray
  """
  if status == "singleton":
    futures = [get_count_table_control.remote(i, subtype, offset, pop, base_dir) for i in range(1,23)]
    results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
  else:
    futures = [get_count_table_control.remote(i, subtype, offset, pop, base_dir, start = 1) for i in range(1,23)]
    results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
  return(results)

def fit_model_all(subtype, offset, pop, base_dir):
  s_tab = get_count_all(subtype, offset, pop, base_dir, status = "singleton")
  s_tab = s_tab.reset_index(level=0)
  s_tab.columns = ['nuc', 'singletons']
  c_tab = get_count_all(subtype, offset, pop, base_dir, status = "control")
  c_tab = c_tab.reset_index(level=0)
  c_tab.columns = ['nuc', 'controls']
  df2 = pd.DataFrame.merge(s_tab, c_tab, on='nuc')
  df = pd.melt(df2, id_vars = 'nuc', var_name = 'status', value_name = "n")
  mod = smf.glm("n ~ status + nuc", df, family = sm.families.Poisson()).fit()
  n_s = sum(s_tab.singletons)
  n_c = sum(c_tab.controls)
  df2['p'] = df2['singletons'] / df2['singletons'].sum()
  df2['q'] = df2['controls'] / df2['controls'].sum()
  df2['tv'] = abs(df2['p'] - df2['q'])
  tot_var = df2['tv'].sum()
  max_d = round(df2['tv'].max(),5)
  mean_var = df2['tv'].mean()
  return {"dev":mod.deviance, "singletons":n_s, "controls":n_c, "offset":offset, "tot_var":tot_var,"mean_var": mean_var ,"max_var": max_d}

ray.init(num_cpus=22)
results = []

print("Running models for subtype: {} in population: {}".format(subtype, pop))
for offset in range(1, 1001):
  print(offset, flush=True)
  results.append(fit_model_all(subtype, offset * -1, pop, base_dir))
  if subtype.startswith("cpg") and offset == 1:
    continue
  results.append(fit_model_all(subtype, offset, pop, base_dir))

ray.shutdown()
  
final = pd.DataFrame.from_dict(results)
out_dir = "{}/single_pos_cc/{}/".format(base_dir, pop)
file_name = out_dir + subtype + "_v2" + ".csv"
final.to_csv(file_name, index = False)

