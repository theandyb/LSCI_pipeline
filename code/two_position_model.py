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
  return pos_list.astype(int)

def s_pos(subtype, chromosome, pop = "ALL"):
  """Get the positions for singletons for a given subtype"""
  input_dir = "output/singletons/{}/pos_files/".format(pop)
  f_name = input_dir + subtype + "_" + str(chromosome) + ".txt"
  pos_list = pd.read_csv(f_name, header=None, names = ['pos'], usecols=['pos']).squeeze("columns")
  return pos_list.astype(int)

def reverse_complement(dna_sequence):
  complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N':'N', 'a':'t', 't':'a', 'c':'g', 'g':'c', 'n':'n'} #added lowercase and N support.
  reverse_seq = dna_sequence[::-1]  # Reverse the sequence
  reverse_complement_seq = "".join(complement.get(base, base) for base in reverse_seq) #get complement, or return original base if not found.
  return reverse_complement_seq

def complement(dna_sequence):
  complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N':'N', 'a':'t', 't':'a', 'c':'g', 'g':'c', 'n':'n'} #added lowercase and N support.
  complement_seq = "".join(complement.get(base, base) for base in dna_sequence) #get complement, or return original base if not found.
  return complement_seq

def generate_nucleotide_motifs(n):
  if n <= 0:
    return []
  nucleotides = ['A', 'C', 'G', 'T']
  motifs = ['']
  for _ in range(n):
    new_motifs = []
    for motif in motifs:
      for nucleotide in nucleotides:
        new_motifs.append(motif + nucleotide)
    motifs = new_motifs
  return {item: 0 for item in motifs}

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
def get_3mer_control(chromosome, pop, subtype, rel_vals=[-1,0,1]):
  motif_size = len(rel_vals)
  control_pos = c_pos(subtype, chromosome, pop)
  rev_control_pos = c_pos(subtype + "_rev", chromosome, pop)
  seq = get_seq_db(chromosome)
  results = generate_nucleotide_motifs(motif_size)
  for index, value in control_pos.items():
    ix_vals = [(value - 1) + item for item in rel_vals]
    nuc = "".join([seq[item] for item in ix_vals])
    if nuc in results.keys():
      results[nuc] += 1
  for index, value in rev_control_pos.items():
    ix_vals = [(value - 1) + (item*-1)  for item in rel_vals]
    ix_vals.sort()
    nuc = "".join([seq[item] for item in ix_vals])
    nuc = reverse_complement(nuc)
    if nuc in results.keys():
      results[nuc] += 1
  return {key: value for key, value in results.items() if value != 0}

@ray.remote
def get_3mer_singletons(chromosome, pop, subtype, rel_vals=[-1,0,1]):
  motif_size = len(rel_vals)
  singleton_pos = s_pos(subtype, chromosome, pop)
  rev_singleton_pos = s_pos(subtype + "_rev", chromosome, pop)
  seq = get_seq_db(chromosome)
  results = generate_nucleotide_motifs(motif_size)
  for index, value in singleton_pos.items():
    ix_vals = [(value - 1) + item for item in rel_vals]
    nuc = "".join([seq[item] for item in ix_vals])
    if nuc in results.keys():
      results[nuc] += 1
  for index, value in rev_singleton_pos.items():
    ix_vals = [(value - 1) + (item*-1)  for item in rel_vals]
    ix_vals.sort()
    nuc = "".join([seq[item] for item in ix_vals])
    nuc = reverse_complement(nuc)
    if nuc in results.keys():
      results[nuc] += 1
  return {key: value for key, value in results.items() if value != 0}

def get_count_all(subtype, status = "singleton", pop = "ALL", rel_vals=[-1,1]):
  """
  Get the singleton and control counts for a given relative position
  across all 22 autosomes. This version is parallelized via ray
  """
  if status == "singleton":
    futures = [get_3mer_singletons.remote(i, pop, subtype, rel_vals=rel_vals) for i in range(1,23)]
    results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
  else:
    futures = [get_3mer_control.remote(i, pop, subtype, rel_vals = rel_vals) for i in range(1,23)]
    results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
  return(results)

def fit_model_all(subtype, p, q, pop = "ALL"):
  df_s = get_count_all(subtype, pop = pop, rel_vals=[p, q])
  df_c = get_count_all(subtype, status = "control", pop = pop, rel_vals=[p, q])
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
  print("p1: {}".format(p1), flush=True)
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
    print("p2: {}".format(p2), flush=True)
    results.append(fit_model_all(subtype, p1, p2, pop = pop))
    final = pd.DataFrame.from_dict(results)
    final.to_csv(file_name, index = False, mode = 'a', header=False)
    results=[]
  start_q = p1+2

ray.shutdown()

