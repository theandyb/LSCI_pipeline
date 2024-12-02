"""
Code for getting the n--mer motif centered at a given list of positions
"""
import pandas as pd
from pyfaidx import Fasta
from Bio.Seq import Seq
import argparse

def get_motif(seq, pos, bp = 10):
    """Get the motif centered at a position from VCF"""
    pyix = pos - 1 # python is index 0, versus 1 for VCF
    return seq[(pyix - bp):(pyix + bp + 1)].seq

def full_cat(ref, alt, motif, bp = 10):
    nuc_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    ref_rc = nuc_dict[ref]
    alt_rc = nuc_dict[alt]
    if ref in ['A','G']:
        final = "{}{}_{}{}".format(ref, ref_rc, alt, alt_rc)
    else:
        final = "{}{}_{}{}".format(ref_rc, ref, alt_rc, alt)
    if motif[bp:bp+2] == "CG":
        final = "cpg_" + final
    return final
  
def full_motif(motif, ref):
    motif_seq = Seq(motif)
    motif_rc = motif_seq.reverse_complement().__str__()
    if ref in ['A','C']:
        final = "{}({})".format(motif, motif_rc)
    else:
        final = "{}({})".format(motif_rc, motif)
    return final

parser = argparse.ArgumentParser(description="Annotate genomic locations with motif")
parser.add_argument("-s", "--singleton", help="Singleton file", required=True)
parser.add_argument("-o", "--output", help="Output file", required=True)
parser.add_argument("-r", "--reference", help="Reference Genome", required=True)
args = parser.parse_args()

singleton_file = args.singleton
output_file = args.output
ref_file = args.reference

print("Reading reference file...")
fasta_obj = Fasta(ref_file)
chrom = "chr1"
current_chrom = chrom
seq = fasta_obj[chrom]

output_list = []

with open(singleton_file) as fp:
    line = fp.readline()
    cnt = 1
    while line:
        data = line.strip().split("\t") # CHROM, POS, REF
        chrom = data[0]
        pos = int(data[1])
        ref = data[2]
        if(chrom != current_chrom):
          current_chrom = chrom
          seq = fasta_obj[chrom]
        motif = get_motif(seq, pos)
        motif_full = full_motif(motif, ref)
        entry = {
            'chrom' : chrom,
            'pos' : pos,
            'ref' : ref,
            'motif_full': motif_full,
        }
        output_list.append(entry)
        cnt += 1
        if cnt % 10000 == 0:
            pd.DataFrame(output_list).to_csv(output_file, index = None, header=False, mode='a')
            output_list = []
        line = fp.readline()

if output_list:
    pd.DataFrame(output_list).to_csv(output_file, index = None, header=False, mode='a')