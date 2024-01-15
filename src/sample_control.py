from pyfaidx import Fasta
import regex as re
import random
import pandas as pd
import argparse
from Bio.Seq import Seq

def main(ref_prefix = "chr"):
  parser = argparse.ArgumentParser(description="Sample control distribution.")
  parser.add_argument("-s", "--singleton", help="Location of singleton file", required=True)
  parser.add_argument("-f", "--fasta", help="FASTA file to grab sequence from", required=True)
  parser.add_argument("-o", "--output", help="Path to output", required=True)
  parser.add_argument("-t", "--subtype", help="Subtype", required = True)
  parser.add_argument("-n", "--nSample", help="Number of controls per singleton", type=int, default = 1)
  args = parser.parse_args()
  ref_file = args.fasta 
  singleton_file = args.singleton 
  nSample = args.nSample
  subtype = args.subtype
  if(subtype[0:3]=="cpg"):
    cpg_bool = True
    all_bool = False
  elif(subtype[0:3]=="all"):
    all_bool = True
    cpg_bool = False
  else:
    cpg_bool = False
    all_bool = False
  current_chrom = 1
  output_list = []
  out_file = args.output
  min_list = []
  max_list = []
  min_file = out_file + ".min"
  max_file = out_file + ".max"
  fasta_obj = Fasta(ref_file)
  seq = fasta_obj["{}{}".format(ref_prefix, current_chrom)]
  seqstr = seq[0:len(seq)].seq
  counter = 1
  with open(singleton_file) as fp:
    line = fp.readline()
    while(line):
      content = line.strip().split("\t") #CHROM, POS, REF
      pos = int(content[1])
      chrom = int(content[0][3:])
      ref = content[2]
      if(current_chrom != chrom):
        current_chrom = chrom
        seq = fasta_obj["{}{}".format(ref_prefix, current_chrom)]
        seqstr = seq[0:len(seq)].seq
      new_entry = sample_control("chr{}".format(current_chrom), pos, ref, nSample, seqstr, cpg_bool = cpg_bool, all_bool = all_bool)
      output_list.extend(new_entry)
      min_dist = min([t.get("distance") for t in new_entry])
      max_dist = max([t.get("distance") for t in new_entry])
      min_entry = [t for t in new_entry if t.get('distance') == min_dist]
      max_entry = [t for t in new_entry if t.get('distance') == max_dist ]
      min_list.extend(min_entry)
      max_list.extend(max_entry)
      line = fp.readline()
      counter += 1
      if counter % 10000 == 0:
        print(counter)
        pd.DataFrame(output_list).to_csv(args.output, index = None, header=False, mode='a')
        output_list = []
        # New: add min and max files
        pd.DataFrame(min_list).to_csv(min_file, index = None, header=False, mode='a')
        min_list = []
        pd.DataFrame(max_list).to_csv(max_file, index = None, header=False, mode='a')
        max_list = []
  print("done sampling")
  if(output_list):
    pd.DataFrame(output_list).to_csv(args.output, index = None, header=False, mode='a')
  if(min_list):
    pd.DataFrame(min_list).to_csv(min_file, index = None, header=False, mode='a')
  if(max_list):
    pd.DataFrame(max_list).to_csv(max_file, index = None, header=False, mode='a')
  return 0

def sample_control(chrom, pos, ref, nSample, seqstr, cpg_bool = False, all_bool = False, window=150, bp=10):
  sites = []
  newlist = []
  off = 0
  if cpg_bool:
    search_str = "CG"
  else:
    if ref == "C" and all_bool == False:
      search_str = "C[ATC]"
    elif ref == "G" and all_bool == False:
      search_str = "[AGT]G"
      off = 1
    else:
      search_str = ref
  while(len(sites) < nSample + 1):
    subseq = seqstr[(pos - 1 - window):(pos+window)]
    sites = [m.start() for m in re.finditer(search_str, subseq, overlapped=True)] # These two lines were changed for CpG/non-CpG sampling
    sites = [s + off for s in sites if (s > bp+window+1 or s < window-bp-1)] # Identify all possible motifs to sample from
    window += 50 #expand window in edge case where mut_site is only ref_allele in window
  window -= 50
  while len(newlist) < nSample:
    if(len(sites)==0):
      print("Bad pos: {}".format(pos))
    ix = random.choice(sites)
    control_al = subseq[ix]
    chrom_ix = ix - window + pos
    #newSeq = seqstr[(chrom_ix - bp - 1):(chrom_ix+bp)].upper()
    #motif2 = full_motif(newSeq, newSeq[bp])
    distance = abs(ix - window)
    entry = {
      'chrom' : chrom,
      'pos' : pos,
      's_ref': ref,
      'c_ref' : control_al,
      'window': window,
      'distance': distance,
      'spos': chrom_ix,
      'three':seqstr[(chrom_ix - 2):(chrom_ix+1)].upper()
    }
    newlist.append(entry)
    sites.remove(ix)
  return newlist

def full_motif(motif, ref):
  motif_seq = Seq(motif)
  motif_rc = motif_seq.reverse_complement().__str__()
  if ref in ['A','C']:
    final = "{}({})".format(motif, motif_rc)
  else:
    final = "{}({})".format(motif_rc, motif)
  return final

if __name__ == "__main__":
  #random.seed( 8675 ) # threeeeee ohhhhh niiiiii-eee-iiiiine
  random.seed( 1776 ) # round two
  main()
