#! /usr/bin/env python
'''
Searches for complexes that could be interesting in transcriptional regulation
using the domain-aware cohesiveness optimization algorithm
'''
import sys, ppi_handler, wDDI_network, gzip
import time

if __name__ == '__main__':
    if len(sys.argv) != 8:
        print 'usage: DACO.py PPIN DDIN TFs PAIR-THRESHOLD PROB-THRESHOLD DEPTH OUT-FILE'
        sys.exit(0)
    
    # ADJUSTMENT: removed data retrieval part and added compatibility for PPIXpress networks
    f_ppin = sys.argv[1]
    f_ddin = sys.argv[2]
    TF_f = sys.argv[3]
    pair_threshold = float(sys.argv[4])
    prob_threshold = float(sys.argv[5])
    depth = int(sys.argv[6])
    out_f = sys.argv[7]
    
    print 'Running domain-aware cohesiveness optimization algorithm using'
    print 'PPIN:', f_ppin
    print 'DDIN:', f_ddin
    print 'TF list:', TF_f
    print 'Pair threshold:', pair_threshold
    print 'Prob. threshold:', prob_threshold
    print 'Max. depth:', depth
    print 'Data-out:', out_f
    print ''
    
    print 'Reading transcription factor file ...'
    tfs = set()
    f_open = open
    if TF_f.endswith("gz") or TF_f.endswith("gzip"):
        f_open = gzip.open
    for line in f_open(TF_f):
        line = line.strip()
        if line != '':
            tfs.add(line)
            
    print 'Reading protein-protein interaction network file ...'
    ppi = ppi_handler.PPI(f_ppin)
    print 'Constructing domain-domain interaction network ... '
    ddi = wDDI_network.DDI_network(f_ddin, ppi)
    
    results = set()
    n = len(tfs)
    start = time.time();
    for i, tf in enumerate(tfs):
        print '============================================================================'
        print str(i+1)+'/'+str(n)+':', tf    # ADJUSTMENT: removed unnecessary output that would have needed name conversion data
        current = ddi.growPairs(tf, tfs = tfs, cutoff = depth, threshold = pair_threshold, prob_threshold = prob_threshold)
        print len(current),'candidates found.'
        sys.stdout.flush()
        results.update(current)
    duration = time.time() - start
    
    # filter proper subsets and complexes that do not contain a seed protein
    to_del = set()
    for i in results:
        if len(i.intersection(tfs)) == 0:
            to_del.add(i)
            continue
        for j in results:
            if i==j: continue
            if i.issubset(j):
                to_del.add(i)
                break
    results.difference_update(to_del)
    
    # write output
    out_f = open(out_f, 'w')
    out_f.writelines([','.join(x)+'\n' for x in results])
    out_f.close()
    
    print "" 
    print len(results),'candidates written to output.'
    print "Overall time:", int(round(duration)), "sec" # added duration