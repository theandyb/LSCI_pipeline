import pandas as pd
import statsmodels.api as sm
#from pyfaidx import Fasta
import sqlite3
from patsy import dmatrices
import statsmodels.formula.api as smf
import argparse
import time
import ray

parser = argparse.ArgumentParser(description="Fit single position models")
parser.add_argument("-p1", "--population", help="1kGP Superpopulation Code", required=True)
parser.add_argument("-p2", "--popComp", help="1kGP Superpopulation Code", required=True)
parser.add_argument("-s", "--subtype", help="Mutation Subtype", required=True)
parser.add_argument("-o", "--output", help="Output base directory", required=True)
parser.add_argument("-c", "--suffix", help="Control file suffix", default="")
args = parser.parse_args()

def c_pos(subtype, chromosome, pop, base_dir, suffix):
  """Get the positions for controls for a given subtype"""
  input_dir = "{}/controls/{}/pos_files/".format(base_dir, pop)
  f_name = input_dir + subtype + "_" + str(chromosome) + ".txt" + suffix
  pos_list = pd.read_csv(f_name, header=None, names = ['pos'], usecols=['pos']).squeeze("columns")
  return pos_list.astype(int)

def s_pos(subtype, chromosome, pop, base_dir):
  """Get the positions for singletons for a given subtype"""
  input_dir = "{}/singletons/{}/pos_files/".format(base_dir, pop)
  f_name = input_dir + subtype + "_" + str(chromosome) + ".txt"
  pos_list = pd.read_csv(f_name, header=None, names = ['pos'], usecols=['pos']).squeeze("columns")
  return pos_list.astype(int)

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
def get_count_table_control(chromosome, subtype, offset, pop, base_dir, suffix):
  """
  Get the count table for controls at a given relative position for a subtype.
  Note: this version has been made so as to be called in parallel via ray, but
  that it sucks in that each worker has to load the reference genome. Need to
  check if there's a python library that has the reference genome as on object
  that can maybe be passed around?
  """
  #fasta_obj = Fasta(ref_file)
  control_pos = c_pos(subtype, chromosome, pop, base_dir, suffix = suffix)
  rev_control_pos = c_pos(subtype + "_rev", chromosome, pop, base_dir, suffix = suffix)
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

@ray.remote
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

def get_count_all(subtype, offset, base_dir, status = "singleton", pop = "ALL", suffix = ""):
  """
  Get the singleton and control counts for a given relative position
  across all 22 autosomes. This version is parallelized via ray
  """
  if status == "singleton":
    futures = [get_count_table_singletons.remote(i, subtype, offset, pop, base_dir) for i in range(1,23)]
    results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
  else:
    futures = [get_count_table_control.remote(i, subtype, offset, pop, base_dir, suffix = suffix) for i in range(1,23)]
    results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
  return(results)

def fit_model(subtype, offset, pop1, pop2, base_dir, suffix=""):
  s_tab1 = get_count_all(subtype, offset, base_dir, status = "singleton", pop = pop1, suffix = suffix)
  s_tab2 = get_count_all(subtype, offset, base_dir, status = "singleton", pop = pop2, suffix = suffix)
  s_tab1 = s_tab1.reset_index(level=0)
  s_tab1.columns = ['nuc', 'singletons']
  s_tab2 = s_tab2.reset_index(level=0)
  s_tab2.columns = ['nuc', 'singletons']
  #s_tab1['pop'] = pop1
  #s_tab2['pop'] = pop2
  
  c_tab1 = get_count_all(subtype, offset, base_dir, status = "control", pop = pop1, suffix = suffix)
  c_tab1 = c_tab1.reset_index(level=0)
  c_tab1.columns = ['nuc', 'controls']
  c_tab2 = get_count_all(subtype, offset, base_dir, status = "control", pop = pop2, suffix = suffix)
  c_tab2 = c_tab2.reset_index(level=0)
  c_tab2.columns = ['nuc', 'controls']
  #c_tab1['pop'] = pop1
  #c_tab2['pop'] = pop2
  
  df1 = pd.DataFrame.merge(s_tab1, c_tab1, on = 'nuc')
  df2 = pd.DataFrame.merge(s_tab2, c_tab2, on = 'nuc')
  df1['pop'] = pop1
  df2['pop'] = pop2
  df = pd.concat([df1, df2])
  df = pd.melt(df, id_vars = ['nuc', 'pop'], var_name = 'status', value_name = "n")
  mod = smf.glm("n ~ status*nuc + nuc*pop + pop*status", df, family = sm.families.Poisson()).fit()
  df['fitted'] = mod.fittedvalues
  df['res'] = mod.resid_deviance
  return(df)

pop1 = args.population
pop2 = args.popComp
subtype = args.subtype
base_dir = args.output
suffix = args.suffix

print("Running with options: ")
print("Populations: " + pop1 + " and " + pop2)
print("Subtype: " + subtype)
print("Base directory: " + base_dir)

time.sleep(3)

ray.init(num_cpus=22)

res_out_dir = "{}/single_pos_pp/{}_{}/".format(base_dir, pop1, pop2)
for offset in range(1, 21):
  print(offset)
  df = fit_model(subtype, offset * -1, pop1, pop2, base_dir, suffix=suffix)
  #fit_model_all(subtype, offset, pop, ref_file, base_dir, suffix = "")
  file_name = res_out_dir + subtype + "_rp_" + str(-1*offset) + ".csv" + suffix
  df.to_csv(file_name, index = False)
  if subtype.startswith("cpg") and offset == 1:
    continue
  df = fit_model(subtype, offset, pop1, pop2, base_dir, suffix=suffix)
  file_name = res_out_dir + subtype + "_rp_" + str(offset) + ".csv" + suffix
  df.to_csv(file_name, index = False)

ray.shutdown()
