"""
Deviance Statistics for Two Position Models
"""
import pandas as pd
import statsmodels.api as sm
import sqlite3
from patsy import dmatrices
import statsmodels.formula.api as smf
import ray
import argparse
import os.path
from sys import exit
#from pyfaidx import Fasta

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--population", help="1kGP Super Population", required=True, choices=["ALL","AFR","AMR","EAS","EUR","SAS"])
parser.add_argument("-t", "--subtype", help="Mutation subtype", choices=["AT_CG", "AT_GC", "AT_TA", "GC_AT", "GC_TA", "GC_CG", "cpg_GC_AT", "cpg_GC_TA", "cpg_GC_CG"], required=True)
parser.add_argument("-s", "--suffix", default = "")
args = parser.parse_args()

pop = args.population
subtype = args.subtype
suffix = args.suffix

def c_pos(subtype, chromosome, pop):
  """Get the positions of control sample locations"""
  input_dir = "output/controls/{}/pos_files/".format(pop)
  f_name = input_dir + subtype + "_" + str(chromosome) + ".txt"
  pos_list = pd.read_csv(f_name, header=0, names = ['pos'], usecols=['pos']).squeeze("columns")
  return pos_list
  
def s_pos(subtype, chromosome, pop = "ALL"):
  """Get the positions for singletons for a given subtype"""
  input_dir = "output/singletons/{}/pos_files/".format(pop)
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

@ray.remote 
def get_count_table_control(chromosome, subtype, p, q, pop):
  """
  Get the count table for nucleotides flanking control locations at two neighboring positions
  """
  control_pos = c_pos(subtype, chromosome, pop)
  rev_control_pos = c_pos(subtype + "_rev", chromosome, pop)
  seqstr = get_seq_db(chromosome)
  # ref_file = "data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
  # fasta_obj = Fasta(ref_file)
  # seq = fasta_obj["{}{}".format("chr", chromosome)]
  # seqstr = seq[0:len(seq)].seq
  results = {"AA":0, "AC":0, "AG":0, "AT":0, 
    "CA":0, "CC":0, "CG":0, "CT":0,
    "GA":0, "GC":0, "GG":0, "GT":0,
    "TA":0, "TC":0, "TG":0, "TT":0}
  for index, value in control_pos.items():
    ix1 = value - 1 + p
    ix2 = value - 1 + q
    nuc = "{}{}".format(seqstr[ix1],seqstr[ix2])
    if nuc in results.keys():
      results[nuc] += 1
  for index, value in rev_control_pos.items():
    ix1 = value - 1 + (p * -1)
    ix2 = value - 1 + (q * -1)
    nuc = "{}{}".format(complement(seqstr[ix1]),complement(seqstr[ix2]))
    if nuc in results.keys():
      results[nuc] += 1
  return(results)

@ray.remote
def get_count_table_control_new(chromosome, subtype, p, q, pop):
  control_pos = c_pos(subtype, chromosome, pop)
  rev_control_pos = c_pos(subtype + "_rev", chromosome, pop)
  seqstr = get_seq_db(chromosome)
  # ref_file = "data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
  # fasta_obj = Fasta(ref_file)
  # seq = fasta_obj["{}{}".format("chr", chromosome)]
  # seqstr = seq[0:len(seq)].seq
  nucleotides_p = [seqstr[pos-1] for pos in (control_pos + p)]
  nucleotides_p_rev = [seqstr[pos-1] for pos in (rev_control_pos + (-1*p))]
  nucleotides_q = [seqstr[pos-1] for pos in (control_pos + q)]
  nucleotides_q_rev = [seqstr[pos-1] for pos in (rev_control_pos + (-1*q))]
  f_all = [a + b for a, b in zip(nucleotides_p, nucleotides_q)]
  rev_all = [a + b for a, b in zip(nucleotides_p_rev, nucleotides_q_rev)]
  frequency_table = {}
  for item in f_all:
    if item in frequency_table:
      frequency_table[item] += 1
    else:
      frequency_table[item] = 1
  for item in rev_all:
    if item in frequency_table:
      frequency_table[item] += 1
    else:
      frequency_table[item] = 1
  return(frequency_table)

@ray.remote
def get_count_table_singletons_new(chromosome, subtype, p, q, pop):
  singleton_pos = s_pos(subtype, chromosome, pop)
  rev_singleton_pos = s_pos(subtype + "_rev", chromosome, pop)
  seqstr = get_seq_db(chromosome)
  # ref_file = "data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
  # fasta_obj = Fasta(ref_file)
  # seq = fasta_obj["{}{}".format("chr", chromosome)]
  # seqstr = seq[0:len(seq)].seq
  nucleotides_p = [seqstr[pos-1] for pos in (singleton_pos + p)]
  nucleotides_p_rev = [seqstr[pos-1] for pos in (rev_singleton_pos + (-1*p))]
  nucleotides_q = [seqstr[pos-1] for pos in (singleton_pos + q)]
  nucleotides_q_rev = [seqstr[pos-1] for pos in (rev_singleton_pos + (-1*q))]
  f_all = [a + b for a, b in zip(nucleotides_p, nucleotides_q)]
  rev_all = [a + b for a, b in zip(nucleotides_p_rev, nucleotides_q_rev)]
  frequency_table = {}
  for item in f_all:
    if item in frequency_table:
      frequency_table[item] += 1
    else:
      frequency_table[item] = 1
  for item in rev_all:
    if item in frequency_table:
      frequency_table[item] += 1
    else:
      frequency_table[item] = 1
  return(frequency_table)

@ray.remote
def get_count_table_singletons(chromosome, subtype, p, q, pop):
  """
  Get the count table for singletons at a given relative position for a subtype.
  """
  singleton_pos = s_pos(subtype, chromosome, pop)
  rev_singleton_pos = s_pos(subtype + "_rev", chromosome, pop)
  seqstr = get_seq_db(chromosome)
  # ref_file = "data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
  # fasta_obj = Fasta(ref_file)
  # seq = fasta_obj["{}{}".format("chr", chromosome)]
  # seqstr = seq[0:len(seq)].seq
  results = {"AA":0, "AC":0, "AG":0, "AT":0, 
    "CA":0, "CC":0, "CG":0, "CT":0,
    "GA":0, "GC":0, "GG":0, "GT":0,
    "TA":0, "TC":0, "TG":0, "TT":0}
  for index, value in singleton_pos.items():
    ix1 = value - 1 + p
    ix2 = value - 1 + q
    nuc = "{}{}".format(seqstr[ix1],seqstr[ix2])
    if nuc in results.keys():
      results[nuc] += 1
  for index, value in rev_singleton_pos.items():
    ix1 = value - 1 + (p * -1)
    ix2 = value - 1 + (q * -1)
    nuc = "{}{}".format(complement(seqstr[ix1]),complement(seqstr[ix2]))
    if nuc in results.keys():
      results[nuc] += 1
  return(results)

# @ray.remote
# def get_count_table_singletons_old(chromosome, subtype, p, q, pop = "ALL"):
#   """
#   Get the count table for singletons at a given relative position for a subtype.
#   Note: this version has been made so as to be called in parallel via ray, but
#   that it sucks in that each worker has to load the reference genome. Need to
#   check if there's a python library that has the reference genome as on object
#   that can maybe be passed around?
#   """
#   singleton_pos = s_pos(subtype, chromosome, pop)
#   rev_singleton_pos = s_pos(subtype + "_rev", chromosome, pop)
#   ref_file = "data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
#   fasta_obj = Fasta(ref_file)
#   seq = fasta_obj["{}{}".format("chr", chromosome)]
#   seqstr = seq[0:len(seq)].seq
#   results = {"AA":0, "AC":0, "AG":0, "AT":0,
#     "CA":0, "CC":0, "CG":0, "CT":0,
#     "GA":0, "GC":0, "GG":0, "GT":0,
#     "TA":0, "TC":0, "TG":0, "TT":0}
#   for index, value in singleton_pos.items():
#     ix1 = value - 1 + p
#     ix2 = value - 1 + q
#     nuc = "{}{}".format(seqstr[ix1],seqstr[ix2])
#     if nuc in results.keys():
#       results[nuc] += 1
#   for index, value in rev_singleton_pos.items():
#     ix1 = value - 1 + (p * -1)
#     ix2 = value - 1 + (q * -1)
#     nuc = "{}{}".format(complement(seqstr[ix1]),complement(seqstr[ix2]))
#     if nuc in results.keys():
#       results[nuc] += 1
#   return(results)

# @ray.remote
# def get_count_table_control_old(chromosome, subtype, p, q, pop = "ALL"):
#   """
#   Get the count table for controls at a pair of given relative positions for a subtype.
#   Note: this version has been made so as to be called in parallel via ray, but
#   that it sucks in that each worker has to load the reference genome. Need to
#   check if there's a python library that has the reference genome as on object
#   that can maybe be passed around?
#   """
#   ref_file = "data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
#   fasta_obj = Fasta(ref_file)
#   control_pos = c_pos(subtype, chromosome, pop)
#   rev_control_pos = c_pos(subtype + "_rev", chromosome, pop)
#   seq = fasta_obj["{}{}".format("chr", chromosome)]
#   seqstr = seq[0:len(seq)].seq
#   results = {"AA":0, "AC":0, "AG":0, "AT":0, 
#     "CA":0, "CC":0, "CG":0, "CT":0,
#     "GA":0, "GC":0, "GG":0, "GT":0,
#     "TA":0, "TC":0, "TG":0, "TT":0}
#   for index, value in control_pos.items():
#     ix1 = value - 1 + p
#     ix2 = value - 1 + q
#     nuc = "{}{}".format(seqstr[ix1],seqstr[ix2])
#     if nuc in results.keys():
#       results[nuc] += 1
#   for index, value in rev_control_pos.items():
#     ix1 = value - 1 + (p * -1)
#     ix2 = value - 1 + (q * -1)
#     nuc = "{}{}".format(complement(seqstr[ix1]),complement(seqstr[ix2]))
#     if nuc in results.keys():
#       results[nuc] += 1
#   return(results)

def get_count_all(subtype, p, q, status = "singleton", pop = "ALL"):
  """
  Get the singleton and control counts for a given relative position
  across all 22 autosomes. This version is parallelized via ray
  """
  if status == "singleton":
    futures = [get_count_table_singletons.remote(i, subtype, p, q, pop) for i in range(1,23)]
    results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
  else:
    futures = [get_count_table_control.remote(i, subtype, p, q, pop) for i in range(1,23)]
    results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
  return(results)

def get_count_all_new(subtype, p, q, status = "singleton", pop = "ALL"):
  """
  Get the singleton and control counts for a given relative position
  across all 22 autosomes. This version is parallelized via ray
  """
  if status == "singleton":
    futures = [get_count_table_singletons_new.remote(i, subtype, p, q, pop) for i in range(1,23)]
    results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
  else:
    futures = [get_count_table_control_new.remote(i, subtype, p, q, pop) for i in range(1,23)]
    results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
  return(results)

# def get_count_all_old(subtype, p, q, status = "singleton", pop = "ALL"):
#   """
#   Get the singleton and control counts for a given relative position
#   across all 22 autosomes. This version is parallelized via ray
#   """
#   if status == "singleton":
#     futures = [get_count_table_singletons_old.remote(i, subtype, p, q, pop) for i in range(1,23)]
#     results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
#   else:
#     futures = [get_count_table_control_old.remote(i, subtype, p, q, pop) for i in range(1,23)]
#     results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
#   return(results)

def fit_model_all(subtype, p, q, pop = "ALL"):
  df_s = get_count_all_new(subtype, p, q, pop = pop)
  df_c = get_count_all_new(subtype, p, q, status = "control", pop = pop)
  df_s = df_s.rename_axis("motif").reset_index()
  df_c = df_c.rename_axis("motif").reset_index()
  df_s.columns = ["motif", "singletons"]
  df_c.columns = ["motif", "controls"]
  df_s = df_s[df_s.motif.str.match("[ACGT]{2}")]
  df_c = df_c[df_c.motif.str.match("[ACGT]{2}")]
  df_s['p1'] = df_s['motif'].apply(lambda x: x[0])
  df_s['p2'] = df_s['motif'].apply(lambda x: x[1])
  df_c['p1'] = df_c['motif'].apply(lambda x: x[0])
  df_c['p2'] = df_c['motif'].apply(lambda x: x[1])
  df_s = df_s.drop(columns = 'motif')
  df_c = df_c.drop(columns = 'motif')
  df = pd.DataFrame.merge(df_s, df_c, on=['p1', 'p2'])
  df = pd.melt(df, id_vars = ['p1', 'p2'], var_name = 'status', value_name = "n")
  mod = smf.glm("n ~ status + p1 + p2 + p1*p2 + p1*status + p2*status", df, family = sm.families.Poisson()).fit()
  n_s = sum(df_s.singletons)
  n_c = sum(df_c.controls)
  dev = mod.deviance.item()
  return {"dev":f'{dev:.2f}', "singletons":n_s, "controls":n_c, "rp1":p, "rp2":q}

print("Running models for subtype: {} in population: {}".format(subtype, pop))
results=[]
out_dir = "output/two_pos/{}/".format(pop)
file_name = out_dir + subtype + ".csv"

if os.path.isfile(file_name):
  with open(file_name) as f:
    for line in f:
      pass
    last_line = line.rstrip().split(",")
  start_p = int(last_line[3])
  start_q = int(last_line[4])
else:
  start_p = -20
  start_q=-19

if start_p == 19 and start_q == 20:
  print("All finished!")
  exit()

ray.init(num_cpus=22)

for p1 in range(start_p, 20):
  print("p1: {}".format(p1))
  if p1 == 0:
    start_q += 1
    continue
  if subtype.startswith("cpg") and p1 == 1:
    start_q += 1
    continue
  for p2 in range(start_q,21):
    if p2 == 0: 
      continue
    if subtype.startswith("cpg") and p2 == 1:
      continue
    print("p2: {}".format(p2))
    results.append(fit_model_all(subtype, p1, p2, pop = pop))
    final = pd.DataFrame.from_dict(results)
    final.to_csv(file_name, index = False, mode = 'a', header=False)
    results=[]
  start_q = p1+2

ray.shutdown()

