#! /usr/bin/env python
'''
Searches for complexes that could be interesting in transcriptional regulation
using the domain-aware cohesiveness optimization algorithm
'''
import sys
from daco_classes import ppi_handler
from daco_classes import wDDI_network
from daco_classes import name_conv

if __name__ == '__main__':
    if len(sys.argv) != 6:
        print 'usage: DACO.py PPIN TFs PAIR-THRESHOLD DEPTH OUT-FILE'
        sys.exit(0)
    f = sys.argv[1]
    TF_f = sys.argv[2]
    pair_threshold = float(sys.argv[3])
    depth = int(sys.argv[4])
    out_f = sys.argv[5]
    
    print 'Running domain-aware cohesiveness optimization algorithm using'
    print 'PPIN:', f
    print 'TF list:', TF_f
    print 'Pair threshold:', pair_threshold
    print 'Max depth:', depth
    print 'Data-out:', out_f
    print ''
    
    print 'Reading transcription factor file ...'
    tfs = set()
    for line in open(TF_f):
        line = line.strip()
        if line != '':
            tfs.add(line)
            
    print 'Reading protein-protein interaction network file ...'
    ppi = ppi_handler.PPI(f)
    print 'Constructing domain-domain interaction network ... '
    ddi = wDDI_network.DDI_network(ppi.getProteins(), ppi)
    
    results = set()
    n = len(tfs)
    for i, tf in enumerate(tfs):
        print '============================================================================'
        print str(i+1)+'/'+str(n)+': Compute seed from', tf,'('+name_conv.NameFromUniprot(tf)+')'
        current = ddi.growPairs(tf, tfs = tfs, cutoff = depth, threshold = pair_threshold)
        print len(current),'candidates found.'
        sys.stdout.flush()
        results.update(current)
    
    # filter proper subsets
    sub = set()
    for i in results:
        for j in results:
            if i==j: continue
            if i.issubset(j):
                sub.add(i)
    results.difference_update(sub)
    
    # write output
    out_f = open(out_f, 'w')
    out_f.writelines([','.join(x)+'\n' for x in results])
    out_f.close()
    print len(results),'candidates written to output.'