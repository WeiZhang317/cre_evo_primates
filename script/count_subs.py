#!/usr/bin/env python

###
# ypar
# in python3.4
# count substitutions
# post count qc includes outlier removal
###

# import modules

import pandas as pd
from sys import argv


def compareseq(exa, exb):
  counts = [x != y if (x is not '-' and x is not '.' and y is not '-' and y is not '.') else '' for (x, y) in zip(exa,exb)].count(True)
  totals = [x == y if (x is not '-' and x is not '.' and y is not '-' and y is not '.') else '' for (x, y) in zip(exa,exb)].count(True)
  return counts, totals


# example

#filename = 'primates/chipseq/peaks/ortho/seqcons/HS/HS_cons_H3K4me1_FDR1_CSpeaks_peak_10258_e39_500bps_seq.txt'

filename = argv[1]

wkdir = 'primates/chipseq/peaks/ortho/seqcons/HS/'
outdir = 'primates/chipseq/peaks/ortho/seqcons/HS_subs/'

outfile = filename.replace('_seq.txt', '_seq_w_subcounts.txt')
outfile = outfile.replace(wkdir, outdir)

peakname = filename.replace(wkdir, '')
peakname = peakname.replace('_seq.txt', '')

sequences = pd.read_csv(filename, sep=' ', header=None, names=['species', 'sequence'])

reference = sequences[sequences.species == 'homo_sapiens']

if len(reference.index) != 1:
  reference = reference.iloc[0]

#comparison = sequences[sequences.species != 'homo_sapiens']

comparison = sequences.groupby('species').first().reset_index()

if not reference.empty and not comparison.empty:
  
  comparison['substitutions'] = 'NA'
  comparison['total_matches'] = 'NA'
  comparison['subs_per_matches'] = 'NA'
  comparison['peakname'] = peakname
  
  for i,r in comparison.iterrows():
    sp = r['species']
    seq = r['sequence']
    
    #print(sp)
    #print(seq)
    
    homs = reference['sequence']
    
    if isinstance(homs, pd.Series):
      homs = homs.values[0]
    else:
      pass
    
    counts, totals = compareseq(homs, seq)
    
    #print(counts)
    #print(totals)
    
    if totals > 0:
      rates = counts / (totals + counts)
    else:
      rates = 'NA'
    
    
    comparison.ix[i, 'substitutions'] = counts
    comparison.ix[i, 'total_matches'] = totals
    comparison.ix[i, 'subs_rate'] = rates
  
  comparison = comparison[['peakname', 'species', 'substitutions', 'total_matches', 'subs_rate']]
  comparison.to_csv(outfile, sep='\t', index=False)


