#! /usr/bin/env python
'''
    Helper class to store PPI-data and provide some little helpers
'''
from collections import defaultdict

class PPI:
    '''
    Helper class to store PPI-data and provide some little helpers
    '''
    def __init__(self, preppi_data, proteins = None, selfIAS = False):
        self._partners = defaultdict(set)
        self._weights = {}
        self._w_whole = defaultdict(float)
        self._selfIAs = selfIAS
        
        if isinstance(preppi_data, str):
            preppi_data = open(preppi_data)
            
        if hasattr(preppi_data, 'read'):
            for line in preppi_data:
                temp = line.strip().split('\t')
                prot1 = temp[0]
                prot2 = temp[1] #added filter for selfIAs
                if (proteins != None and prot1 in proteins and prot2 in proteins and (selfIAS or prot1 != prot2)) or proteins == None:
                    weight = float(temp[2])
                    self._partners[prot1].add(prot2)
                    self._partners[prot2].add(prot1)
                    self._weights[(prot1, prot2)] = weight
                    self._weights[(prot2, prot1)] = weight
                    self._w_whole[prot1] += weight
                    self._w_whole[prot2] += weight
        else: #preppi_data is PPI instance
            for prot1 in preppi_data._partners:
                if proteins != None and prot1 not in proteins:
                    continue
                for prot2 in preppi_data._partners[prot1]:
                    if (proteins != None and prot2 not in proteins) or prot1 < prot2: #only included once and no double-enumeration
                        continue
                    self._partners[prot1].add(prot2)
                    self._partners[prot2].add(prot1)
                    weight = preppi_data._weights[(prot1, prot2)]
                    self._weights[(prot1, prot2)] = weight
                    self._weights[(prot2, prot1)] = weight
                    
    def getSubnetwork(self, proteins):
        '''
        Returns subnetwork of PPI-instance only containing data about the given proteins
        '''
        return PPI(self, proteins)
    
    def size(self):
        '''
        Number of proteins within the network
        '''
        return len(self._partners.keys())
    
    def writeToFile(self, filename, selfIA = True):
        '''
        Writes tab-separated network to given path
        '''
        lines = ['Protein1\tProtein2\tweight\n']
        
        for prot1 in self._partners:
            for prot2 in self._partners[prot1]:
                if prot1 < prot2:
                    continue
                if prot1 == prot2 and selfIA == False:
                    continue
                lines.append('\t'.join([prot1, prot2, str(self._weights[(prot1, prot2)])+'\n']))
                
        open(filename, 'w').writelines(lines)
    
    def computeCohesiveness(self, proteins):
        '''
        Computes the cohesiveness for a set of proteins (UniProt-IDs) within the given PPI-network,
        weights_within/(weights_within+weights_to_outside)
        '''
        if len(proteins) <= 1:
            return 0.0
        w_in = 0.0
        w_out = 0.0
        for prot1 in proteins: # do not allow for homo-oligomers
            for prot2 in (x for x in self._partners[prot1] if x != prot1):
                # every edge is only counted once 
                if prot2 in proteins and prot1 < prot2:
                    w_in += self._weights[(prot1, prot2)]
                else: # is naturally counted once since we enumerate "from the inside"
                    w_out += self._weights[(prot1, prot2)]
                    
        return w_in / (w_out + w_in)
    
    def computeClusterCohesiveness(self, proteins):
        '''
        Computes the cohesiveness for a set of proteins (UniProt-IDs) within the given PPI-network, but outputs components
        (weights_within, weights_to_outside)
        '''
        w_in = 0.0
        w_out = 0.0
        for prot1 in proteins: # do not allow for homo-oligomers
            for prot2 in (x for x in self._partners[prot1] if x != prot1):
                # every edge is only counted once 
                if prot2 in proteins and prot1 < prot2:
                    w_in += self._weights[(prot1, prot2)]
                else: # is naturally counted once since we enumerate "from the inside"
                    w_out += self._weights[(prot1, prot2)]
                    
        return (w_in, w_out)
    
    def computeDeltaCohesiveness(self, protein, proteins):
        '''
        Computes the cohesiveness for a proteins (UniProt-IDs) within the given PPI-network and component proteins,
        output (weights_within, weights_to_outside) for this protein
        '''
        w_in = 0.0
        for prot2 in proteins:
            if prot2 in self._partners[protein]:
                w_in += self._weights[(protein, prot2)]
                    
        return (w_in, self._w_whole[protein] - w_in)
    
    def getProteins(self):
        '''
        Returns the set of proteins involved in the protein-protein interaction network
        '''
        return (set(self._partners.keys()))
    
    def getProb(self, proteins):
        '''
        Computes the probability of the chain p1, ..., pn
        '''
        P = 1.0
        for i in range(len(proteins)-1):
            P *= self._weights(proteins[i],proteins[i+1])
        
        return P