from pyfaidx import Fasta
import argparse
import pandas as pd
from functools import partial, partialmethod

pd.DataFrame.to_tsv = partialmethod(pd.DataFrame.to_csv, sep='\t', index = None)

parser = argparse.ArgumentParser(description="Identify CpG and Non-CpG singletons")
parser.add_argument("-s", "--singleton", help="Location of singleton file", required=True)
parser.add_argument("-o", "--output", help="Location of output file", required=True)
parser.add_argument("-r", "--reference", help="Reference genome fasta", required=True)
args = parser.parse_args()

singleton_file = args.singleton
output_file = args.output
ref_file = args.reference

fasta_obj = Fasta(ref_file)

current_chrom = 1
seq = fasta_obj["chr{}".format(current_chrom)]
seqstr = seq[0:len(seq)].seq

output_list = []

with open(singleton_file) as fp:
    cnt = 1
    line = fp.readline()
    while(line):
        data = line.strip().split("\t")
        chrom = int(data[0][3:])
        pos = int(data[1])
        ref_al = data[2]
        if(current_chrom != chrom):
            current_chrom = chrom
            seq = fasta_obj["chr{}".format(current_chrom)]
            seqstr = seq[0:len(seq)].seq
        ref = seqstr[(pos-1):(pos+1)]
        ref2 = seqstr[(pos-2):(pos)]
        if(ref == "CG" or ref2 == "CG"):
            cpg_stat = 1
        else:
            cpg_stat = 0
        entry = {'chrom': chrom, 'pos': pos, 'ref': ref_al, 'cpg': cpg_stat}
        output_list.append(entry)
        cnt += 1
        if cnt % 100000 == 0:
            print(cnt)
            pd.DataFrame(output_list).to_tsv(output_file, header = False, mode='a')
            output_list = []
        line = fp.readline()

if output_list:
    pd.DataFrame(output_list).to_tsv(output_file, header = False, mode='a')