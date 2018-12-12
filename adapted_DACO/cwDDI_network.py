#! /usr/bin/env python
import sys, gzip
from collections import defaultdict
'''
weighted DDI-network class
'''

def PSort(t):
        '''
        Sorts pair fast and returns [sorted-pair]
        '''
        if t[0] < t[1]:
            return [t]
        else:
            return [(t[1], t[0])]
        
class Protein:
    '''
    Stores one single protein, its gene name, Pfam-predicted domain(s)
    '''
    def __init__(self, uniprot_identifier): # added reading from file only, completely stripped retrieval
        self._up_id = uniprot_identifier
        
        # get possible _domains
        self._domains = []
        
    def init_domain_multiplicity(self):
        '''
        Look for multiplicity of a domain and save mapping to domain instances in self._mult_pfam_map
        '''
        pfam_count = defaultdict(int)
        self._mult_pfam_map = defaultdict(list)
        for dom in self._domains:
            pfam_count[dom._pfam_id] += 1
            self._mult_pfam_map[dom._pfam_id].append(dom)
        for pfam_id in pfam_count:
            if pfam_count[pfam_id] == 1:
                del self._mult_pfam_map[pfam_id]
    
    def isDisconnected(self):
        '''
        Determines if node is disconnected from the remaining network
        '''
        for domain in self._domains:
            if domain.hasInteraction():
                return False
        return True
    
    def cleanup(self):
        '''
        Removes domains not involved in any interaction
        '''
        
        fb = None
        reachable = set()
        for d in self._domains:
            if d._pfam_id == 'FB':
                fb = d
                continue
            else:
                for d2 in d._partners:
                    reachable.add(d2._protein)
        
        # delete directly reachable ones from fb if fb are used
        if fb != None:
            for fb2 in fb._partners[:]:
                if fb2._protein in reachable:
                    fb2._partners.remove(fb)
                    fb._partners.remove(fb2)
                    
        self._domains = [x for x in self._domains if x.hasInteraction()]
        
    def __str__(self):
        return self._up_id+'|'+self._gene_name+'('+','.join((str(x) for x in self._domains))+')'
    
    def __repr__(self):
        return self._up_id
    
    def saturatedDomains(self):
        return sum([1 for x in self._domains if x._isSaturated])
    
    def __hash__(self):
        return self._up_id.__hash__()
    
    def __eq__(self, other):
        return self._up_id.__eq__(other)
    
    def __ne__(self, other):
        return self._up_id.__ne__(other)
    
    def __gt__(self, other):
        return self._up_id.__gt__(other)
    
    def __le__(self, other):
        return self._up_id.__le__(other)
    
    def __lt__(self, other):
        return self._up_id.__lt__(other)
    
class Domain:
    '''
    Stores domain and domain-interaction information
    '''
    def __init__(self, pfam_identifier, protein, number):
        self._pfam_id = pfam_identifier
        self._protein = protein
        self._domain_number = number
        self._partners = []
        self._name = '|'.join([str(self._domain_number), self._pfam_id, self._protein._up_id])
        protein._domains.append(self)
        
    def hasInteraction(self):
        return len(self._partners) > 0
    
    def __str__(self):
        return self._name
    
    def __repr__(self):
        return self._name
    
    def __hash__(self):
        return self._name.__hash__()
    
    def __eq__(self, other):
        return self._name.__eq__(other)
    
    def __ne__(self, other):
        return self._name.__ne__(other)
    
    def __gt__(self, other):
        return self._name.__gt__(other)
    
    def __le__(self, other):
        return self._name.__le__(other)
    
    def __lt__(self, other):
        return self._name.__lt__(other)
         
class DDI_network:
    '''
    Builds and stores DDI-network given a list of UniProt-ids and a PPI-object, 
    DDI are only established if there is an edge in the PPIN with a weight > threshold.
    '''
    def __init__(self, ddi_path, ppi): # added reading from file only, completely stripped retrieval
        self._proteins = []
        self._ppi = ppi
        self._up_map = {}
        self._amb_remove = 0
        self._edge_grow_hash = set()
        self._prob_threshold = 0.5 # added to accomodate
        
        f_open = open
        if ddi_path.endswith(".gz") or ddi_path.endswith(".gzip"):
            f_open = gzip.open
        
        # ADJUSTMENT: PPIXpress compatibility
        domain_map = {}
        for line in f_open(ddi_path):
            if line.startswith("Protein/Domain1"):
                continue
            t1, t2, IA_type, _ = line.strip().split()
            if IA_type == 'pd': # protein-domain
                if not t1 in self._up_map:
                    p = Protein(t1)
                    self._up_map[p._up_id] = p
                else:
                    p = self._up_map[t1]
                d_spl = t2.split("|")
                d = Domain(d_spl[1], p, int(d_spl[0]))
                domain_map[t2] = d
            else: # domain-domain
                d1 = domain_map[t1]
                d2 = domain_map[t2]
                d1._partners.append(d2)
                d2._partners.append(d1)
            
        del domain_map
        
        # call multiplicity detection
        for p in self._proteins:
            p.init_domain_multiplicity()
            
        # remove disconnected proteins from network
        self.cleanup()
        for prot in self._proteins:
            self._up_map[prot._up_id] = prot
            
    def cleanup(self):
        '''
        Removes disconnected proteins and unneeded domains
        '''
        self._proteins = [x for x in self._proteins if not x.isDisconnected()]
        for protein in self._proteins: protein.cleanup()
        
    def __len__(self):
        return len(self._proteins)
                                                   
    def __str__(self):
        return 'DDI of '+str(self._proteins)
    
    def __repr__(self):
        return 'DDI of '+str(self._proteins)
    
    def writeDDIN(self, path, prots_of_interest = []):
        '''
        Writes a tab-separated file with format PROTEIN1/DOMAIN1 PROTEIN2/DOMAIN2 ia type, where weights between _proteins and their _domains are very high
        '''
        lines = []
        lines.append('\t'.join(['PROT1/DOM1', 'PROT2/DOM2', 'weight', 'type']))
        for prot in self._proteins:
            for domain in prot._domains:
                if len(prots_of_interest) == 0 or prot._up_id not in prots_of_interest:
                    lines.append('\t'.join([prot._up_id, str(domain), '1', 'pd']))
                else:
                    lines.append('\t'.join([prot._up_id, str(domain), '1', 'pd']))
                for partner in domain._partners:
                    if partner._protein._up_id < domain._protein._up_id:
                        lines.append('\t'.join([str(domain), str(partner),'0.005', 'dd']))
        
        open(path,'w').writelines('\n'.join(lines))
    
    def writeDDIguidedPPIN(self, path, prots_of_interest = []):
        '''
        Writes a tab-separated file with format PROTEIN1 PROTEIN2 #DDI-WEIGHT
        '''
        lines = []
        lines.append('\t'.join(['PROTEIN1', 'PROTEIN2', '#INTERACTIONS']))
        for i in range(len(self._proteins)-1):
            prot1 = self._proteins[i]
            for j in range(i+1, len(self._proteins)):
                prot2 = self._proteins[j]
                ias = 0
                for domain in prot1._domains:
                    for partner in domain._partners:
                        if partner._protein == prot2:
                            ias +=1
                if ias > 0:
                    if len(prots_of_interest)==0 or (prot1._up_id not in prots_of_interest and prot2._up_id not in prots_of_interest):
                        lines.append('\t'.join([prot1._up_id, prot2._up_id, str(ias)]))
                    else:
                        lines.append('\t'.join([prot1._up_id, prot2._up_id, str(ias)]))
                
                        
        open(path,'w').writelines('\n'.join(lines))  
    
    def writePPISubnetwork(self, out_path):
        '''
        Writes the PrePPI-network of the involved proteins to out_path:
        tab-separated: PROTEIN1 PROTEIN2 PrePPI-weight
        '''
        PPI_subnetwork =self._ppi.getSubnetwork(set([x._up_id for x in self._proteins]))
        PPI_subnetwork.writeToFile(out_path)
    
    def getIncidentNodes(self, proteins, occupied_domains):
        '''
        Given a set/list of proteins (UPId) and occupied domains, return map DDI-reachable proteins->domain-pair used
        '''
        adjacent = defaultdict(set)
        for prot in (self._up_map[x] for x in proteins):
            for d in (x for x in prot._domains if x not in occupied_domains):
                for d2 in (x for x in d._partners if (x not in occupied_domains and x._protein._up_id not in proteins)):
                    adjacent[d2._protein._up_id].add((d, d2))
                    
        return adjacent
    
    def SafeGetIncidentNodes(self, proteins, occupied_domains = set()):
        '''
        Given a set/list of proteins (UPId) and occupied domains, return map DDI-reachable proteins->domain-pair used;
        additionally checks for those annotated
        '''
        adjacent = defaultdict(set)
        for prot in (self._up_map[x] for x in proteins if x in self._up_map):
            for d in (x for x in prot._domains if x not in occupied_domains):
                for d2 in (x for x in d._partners if (x not in occupied_domains and x._protein._up_id not in proteins)):
                    adjacent[d2._protein._up_id].add((d, d2))
                    
        return adjacent
    
    def getUsableBoundaryNodes(self, domain_edges):
        '''
        Find usable boundary (incident to border) nodes in proteins
        that have exactly one occ. domain.
        If they had more, graph could be disconnected internally and invariants don't hold anymore.
        '''
        domains_used_map = defaultdict(list)
        dd_m = {}
        for d1, d2 in domain_edges:
            domains_used_map[d1._protein._up_id].append(d1)
            domains_used_map[d2._protein._up_id].append(d2)
            dd_m[d1] = d2
            dd_m[d2] = d1
            
        associated_domains = {}
        for prot in (prot for prot in domains_used_map if len(domains_used_map[prot]) == 1):
            domain = domains_used_map[prot][0]
            associated_domains[prot] = (domain, dd_m[domain])
                
        return associated_domains
    
    def filterAlternatives(self, pairs):
        '''
        Filters tupel of domain-choices to tame combinatorial explosion.
        Pairs is a set of (domain_inside , domain_outside).
        '''
        already_seen = set()
        filtered = []
        for a1, a2 in pairs:
            hash_str = a1._pfam_id + a1._protein._up_id + a2._pfam_id
            if hash_str in already_seen:
                continue
            
            already_seen.add(hash_str)
            filtered.append( (a1, a2) )
            
        return filtered
        
    def growClusterStep(self, proteins, domain_edges, c_w_in, c_w_out, P):
        '''
        DDI-based clustering implementation
        '''
        c_sum = c_w_in + c_w_out
        max_coh = c_w_in / c_sum
        
        n_w_in = 0.0
        n_w_out = 0.0
        
        # determine occupied domains from DDI-edges
        occupied_domains = set()
        for d1, d2 in domain_edges:
            occupied_domains.add(d1)
            occupied_domains.add(d2)
            
        # search for an add. protein that max cohesiveness
        boundary_new = self.getIncidentNodes(proteins, occupied_domains)
        add_max_prot = None
        for prot in boundary_new:
            (w_in, w_out) = self._ppi.computeDeltaCohesiveness(prot, proteins)
            temp_coh = (c_w_in + w_in) / (c_sum + w_out)
            if temp_coh > max_coh:
                max_coh = temp_coh
                add_max_prot = prot
                n_w_in = w_in
                n_w_out = w_out
                
        # search for removable protein that maximizes current_coh
        boundary_old = self.getUsableBoundaryNodes(domain_edges)
        del_max_prot = None
        for prot in boundary_old:
            (w_in, w_out) = self._ppi.computeDeltaCohesiveness(prot, proteins)
            temp_coh = (c_w_in - w_in) / (c_sum - w_out)
            if temp_coh > max_coh:
                max_coh = temp_coh
                del_max_prot = prot
                n_w_in = w_in
                n_w_out = w_out
        
        # deletion couldnt maximize further -> add or nothing
        if del_max_prot == None:
            # add couldnt max further -> done
            if add_max_prot == None:
                return frozenset(proteins)
            else:
                
                # filter alternatives
                less_alternatives = self.filterAlternatives(boundary_new[add_max_prot])
                if len(less_alternatives) == 1:
                    P_n = P * self._ppi._weights[(less_alternatives[0][0]._protein._up_id, less_alternatives[0][1]._protein._up_id)]
                    if P_n >= self._prob_threshold:
                        return (proteins.union([add_max_prot]), domain_edges.union(PSort(less_alternatives[0])), c_w_in + n_w_in, c_w_out + n_w_out - n_w_in, P_n)
                    else:
                        return frozenset(proteins)
                else:
                    results = []
                    prot_new = proteins.union([add_max_prot])
                    d_in = c_w_in + n_w_in
                    d_out = c_w_out + n_w_out - n_w_in
                    for d in less_alternatives:
                        P_n = P * self._ppi._weights[(d[0]._protein._up_id, d[1]._protein._up_id)]
                        if P_n >= self._prob_threshold:
                            results.append((prot_new, domain_edges.union(PSort(d)), d_in, d_out, P_n))
                    if len(results) == 0:
                        return frozenset(proteins)
                       
                    return results
        else:
            # remove protein and delete related domain-occupation
            d1, d2 = boundary_old[del_max_prot]
            P_n = P / self._ppi._weights[(d1._protein._up_id, d2._protein._up_id)]
            return (proteins.difference([del_max_prot]), domain_edges.difference(PSort(boundary_old[del_max_prot])), c_w_in - n_w_in, c_w_out - n_w_out + n_w_in, P_n)
        
    def growManager(self, cutoff, job):
        '''
        Manages the iterative search, in weighted DDIs variant no fast-mode possible
        '''
        out = set()
        stack = [job]
        
        while len(stack) > 0:
            proteins, domain_edges, c_w_in, c_w_out, P = stack.pop()
            
            str_hash = frozenset(domain_edges)
            
            # already computed -> no duplicate computation needed
            if str_hash in self._edge_grow_hash:
                continue

            # catch cutoff here
            if len(proteins) >= cutoff:
                out.add(frozenset(proteins)) # don't add hash, minimal slower, much less memory?
            else:
                result = self.growClusterStep(proteins, domain_edges, c_w_in, c_w_out, P)
                self._edge_grow_hash.add(str_hash)

                # output is a result
                if isinstance(result, frozenset):
                    out.add(result)
                # a list of alternative steps
                elif isinstance(result, list):
                    stack.extend(result)
                # or one alternative step
                else:
                    stack.append(result)
                    
        return out
        
    def combineGrowthOutput(self, out):
        '''
        Merges a tree-like list to a duplicate-free list of clusters (frozensets)
        '''
        res = set()
        for i in out:
            if isinstance(i, list):
                for i in self.combineGrowthOutput(i):
                    res.add(i)
            else:
                res.add(i)
                
        return list(res)
    
    # ADJUSTMENT: added prob_threshold (was fixed at 0.5)
    def growPairs(self, protein, tfs, cutoff, threshold = 0.9, prob_threshold = 0.5, verbose = False):
        '''
        First builds candidate starting seeds, then seeds.
        DDI-information for starting clusters is partially neglected by allowing connections without domains!
        '''
        w = self._ppi._weights
        self._prob_threshold = prob_threshold # added to accomodate
        
        # needed to determine if in DDIN
        ref_ddi = self.SafeGetIncidentNodes([protein])
        
        if protein in tfs:
            to_del = []
            for p2 in ref_ddi:
                if p2 in tfs and p2 > protein: # ADJUSTMENT: reversed the direction here to have exactly the same input set as in JDACO
                    to_del.append(p2)
            for p2 in to_del:
                del ref_ddi[p2]
            del to_del
            
        tfs = None
        
        # determine max/above threshold
        # all ddi partners
        ddi_partners = ref_ddi.keys()
        
        # pre_sort all DDII-partners to gather some information
        ddi_partners.sort(key=lambda x:w[(protein, x)], reverse = True)
        
        # check if in PPIN
        if len(ddi_partners) == 0:
            'Protein not in PPIN, no seeding possible.'
            return set()
        
        # ADJUSTMENT: differs from original implementation here, that would have taken two partners regardless of threshold
        ddi_partners = [x for x in ddi_partners if w[(x, protein)] > threshold]
        no_pairs = len(ddi_partners)
        
        with_ddi = {}
        for ddi_partner in ddi_partners:
            with_ddi[ddi_partner] = ref_ddi[ddi_partner]
        
        del ddi_partners
        
        alternatives = [(x, self.filterAlternatives(with_ddi[x])) for x in with_ddi]
        del with_ddi
           
        jobs = []
        for partner, variants in alternatives:
            (w_in, w_out) = self._ppi.computeClusterCohesiveness([protein, partner])
            jobs.extend([(set([protein, partner]),set(PSort(d)), w_in, w_out, self._ppi._weights[(protein,partner)]) for d in variants])
        
        results = set()
        no_jobs = len(jobs)
        print 'Seeding from', no_pairs,'pair(s) with', no_jobs, 'unique domain interactions.'
        for i, job in enumerate(jobs):
            result = self.growManager(cutoff, job)
            results.update(result)
            if verbose:
                print 'Finished',str(i+1)+'/'+str(no_jobs) +'.'
                print result
                sys.stdout.flush()
            
        return results