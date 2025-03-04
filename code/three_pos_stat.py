"""
Deviance statistics for three position models
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

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--population", help="1kGP Super Population", required=True, choices=["ALL","AFR","AMR","EAS","EUR","SAS"])
parser.add_argument("-t", "--subtype", help="Mutation subtype", choices=["AT_CG", "AT_GC", "AT_TA", "GC_AT", "GC_TA", "GC_CG", "cpg_GC_AT", "cpg_GC_TA", "cpg_GC_CG"], required=True)
parser.add_argument("-o", "--output", help="Output base directory", required=True)
parser.add_argument("-s", "--suffix", default = "")
args = parser.parse_args()

pop = args.population
subtype = args.subtype
suffix = args.suffix
base_dir = args.output

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

def reverse_complement(dna_sequence):
  complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N':'N', 'a':'t', 't':'a', 'c':'g', 'g':'c', 'n':'n'} #added lowercase and N support.
  reverse_seq = dna_sequence[::-1]  # Reverse the sequence
  reverse_complement_seq = "".join(complement.get(base, base) for base in reverse_seq) #get complement, or return original base if not found.
  return reverse_complement_seq

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
def get_count_table_control(chromosome, subtype, rel_vals = [-1,0,1], pop = "ALL"):
  control_pos = c_pos(subtype, chromosome, pop)
  rev_control_pos = c_pos(subtype + "_rev", chromosome, pop)
  seq = get_seq_db(chromosome)
  results = generate_nucleotide_motifs(3)
  for index, value in control_pos.items():
    ix_vals = [(value - 1) + item for item in rel_vals]
    nuc = "".join([seq[item] for item in ix_vals])
    if nuc in results.keys():
      results[nuc] += 1
  for index, value in rev_control_pos.items(): # reverse complement: need to take inverse (x-1) of each offset, sort, get nucleotide motif, and then reverse-complement
    ix_vals = [(value - 1) + (item*-1)  for item in rel_vals]
    ix_vals.sort()
    nuc = "".join([seq[item] for item in ix_vals])
    nuc = reverse_complement(nuc)
    if nuc in results.keys():
      results[nuc] += 1
  return {key: value for key, value in results.items() if value != 0}

@ray.remote
def get_count_table_singletons(chromosome, subtype, rel_vals = [-1, 0, 1], pop = "ALL"):
  """
  Get the count table for singletons at a given relative position for a subtype.
  Note: this version has been made so as to be called in parallel via ray, but
  that it sucks in that each worker has to load the reference genome. Need to
  check if there's a python library that has the reference genome as on object
  that can maybe be passed around?
  """
  singleton_pos = s_pos(subtype, chromosome, pop)
  rev_singleton_pos = s_pos(subtype + "_rev", chromosome, pop)
  seq = get_seq_db(chromosome)
  results = generate_nucleotide_motifs(3)
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
  return(results)

def get_count_all(subtype, rel_vals = [-1,0,1], status = "singleton", pop = "ALL"):
  """
  Get the singleton and control counts for a given relative position
  across all 22 autosomes. This version is parallelized via ray
  """
  if status == "singleton":
    futures = [get_count_table_singletons.remote(i, subtype, rel_vals = rel_vals, pop = pop) for i in range(1,23)]
    results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
  else:
    futures = [get_count_table_control.remote(i, subtype, rel_vals = rel_vals, pop = pop) for i in range(1,23)]
    results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
  return(results)

def fit_model_all(subtype, rel_vals, pop = "ALL"):
  df_s = get_count_all(subtype, rel_vals = rel_vals, pop = pop)
  df_c = get_count_all(subtype, rel_vals = rel_vals, status = "control", pop = pop)
  df_s = df_s.rename_axis("motif").reset_index()
  df_c = df_c.rename_axis("motif").reset_index()
  df_s.columns = ["motif", "singletons"]
  df_c.columns = ["motif", "controls"]
  df_s['p1'] = df_s['motif'].apply(lambda x: x[0])
  df_s['p2'] = df_s['motif'].apply(lambda x: x[1])
  df_s['p3'] = df_s['motif'].apply(lambda x: x[2])
  df_c['p1'] = df_c['motif'].apply(lambda x: x[0])
  df_c['p2'] = df_c['motif'].apply(lambda x: x[1])
  df_c['p3'] = df_c['motif'].apply(lambda x: x[2])
  df_s = df_s.drop(columns = 'motif')
  df_c = df_c.drop(columns = 'motif')
  df = pd.DataFrame.merge(df_s, df_c, on=['p1', 'p2', 'p3'])
  df = pd.melt(df, id_vars = ['p1', 'p2', 'p3'], var_name = 'status', value_name = "n")
  mod = smf.glm("n ~ status + p1 + p2 + p1*p2*p3 + p1*p2*status + p2*p3*status + p1*p3*status", df, family = sm.families.Poisson()).fit()
  n_s = sum(df_s.singletons)
  n_c = sum(df_c.controls)
  return {"dev":mod.deviance, "singletons":n_s, "controls":n_c, "rp1":rel_vals[0], "rp2":rel_vals[1], "rp3": rel_vals[2]}

ray.init(num_cpus=22)
results = []
print("Running models for subtype: {} in population: {}".format(subtype, pop))

for p1 in range(-10,9):
    if p1 == 0:
        continue
    if subtype.startswith("cpg") and p1 == 1:
        continue
    print("p1: {}".format(p1), flush=True)
    for p2 in range((p1+1),10):
        if p2 == 0:
            continue
        if subtype.startswith("cpg") and p2 == 1:
            continue
        print("p2: {}".format(p2), flush=True)
        for p3 in range((p2+1),11):
            if p3 == 0:
                continue
            if subtype.startswith("cpg") and p3 == 1:
                continue
            print("p3: {}".format(p3), flush=True)
            results.append(fit_model_all(subtype, [p1, p2, p3], pop = pop))

ray.shutdown()
final = pd.DataFrame.from_dict(results)
file_name = base_dir + "/three_pos/" + pop + "/" + subtype + ".csv"
final.to_csv(file_name, index = False)
