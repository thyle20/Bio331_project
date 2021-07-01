# Bio331_project
Alignment of HPV strains

#import matplotlib.pyplot as plt
#from itertools import combinations
import re
from statistics import mode

def read_fasta (file):
  F = open(file)
  seq_dict = {}
  id = ""

  for line in F:
    line = line.rstrip("\n")
    x = re.search(r">(\S+)",line)
    if x:
      id = x.group(1)
      seq_dict[id] = ""
    else:
      seq_dict[id] += line

  F.close()
  return (seq_dict)

gencode = {
  'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
  'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
  'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
  'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
  'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
  'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
  'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
  'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
  'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
  'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
  'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
  'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
  'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
  'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
  'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
  'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

def diversity_output (list):
  hetero_freq = []
  freq_dict = {}
  for value in list:
    freq_dict[value] = round(list.count(value)/len(list),2)
  print (freq_dict, end = "\t")
  unique_value = []
  for key in freq_dict:
    unique_value.append(key)
  print (unique_value, end ="\t")
  hetero = 0
  for index in range (0, len(unique_value)-1):
    val = unique_value[index]
    for index1 in range (index+1, len(unique_value)):
      val1 = unique_value[index1]
      hetero += 2* freq_dict[val]* freq_dict[val1]
  print (round(hetero,2), end = "\t")
  hetero_freq.append(hetero)


file = "HPV_Align.fas"
seq_dict = read_fasta (file)
seq_list = []
for key in seq_dict:
  seq = seq_dict[key]
  seq_list.append(seq)

aa_for_graph =[]
aa_freq_list = []

nuc_for_graph = []
nuc_freq_list = []

cod_variant = []
tmp_codon = []
aa_hetero_freq =[]
cod_hetero_freq = []

i = 0
while i < min(len(index) for index in seq_list):
  print(i,end="\t")
  codon_list = []
  for seq in seq_list:
    codon = seq[i:i+3]
    codon_list.append(codon)

  gap_count = 0

  for x in range(0,3):
    for codon in codon_list:
      nuc = codon[x:x+1]
      if nuc == "-":
        gap_count += 1
    
  if gap_count ==0:
    aa_list = []
    for cod in codon_list:
      amino_acid = gencode[cod]
      aa_list.append(amino_acid)
  
    aa_for_graph.append(mode(aa_list))
    aa_freq_list.append(aa_list.count(mode(aa_list)))

    nuc_for_graph.append(mode(codon_list))
    nuc_freq_list.append(codon_list.count(mode(codon_list)))

    diversity_output(codon_list)
    diversity_output(aa_list)
    
  i += 3
  print()

pos = []
for position in range(0,len(aa_for_graph)):
  pos.append(position)
#print(pos)
#print(nuc_for_graph)
#print(aa_for_graph)
#print(hetero_freq)

#Code for drawing the graph
#plt.plot(pos, nuc_freq_list, label = "codon" )
#plt.plot(pos, aa_freq_list, label = "aa")
#plt.locator_params(axis= "x", nbins = 11)
#plt.xlabel("Position in the sequence")
#plt.ylabel("Frequency")
#plt.title("Graph")
#plt.legend(loc = "right",bbox_to_anchor = (1.14, 0.5), ncol =1, fontsize = 8.5)
#plt.show()
