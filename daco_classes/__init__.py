#! /usr/bin/env python
import name_conv, domain_query, uniprot_handler, ppi_handler, sys
from collections import defaultdict
'''
DDI-network builder and MEI-finder
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
    def __init__(self, uniprot_identifier, use_fallback_domains = True):
        '''
        Builds domain-info from UP-ID, uses Pfam and InterPro.
        Pfam _domains are primary source (multiplicity!), InterPro-Pfam associations are added only if domain not in Pfam
        '''
        self._up_id = uniprot_identifier
        self._gene_name = name_conv.NameFromUniprot(uniprot_identifier)
        
        # get possible _domains
        self._domains = []
        # Pfam domain from Pfam database (NOT uniprot!)
        pfam_domains = [x[0] for x in domain_query.getPfamData(uniprot_identifier)]
        # InterPro annotation from UniProt
        up_ipro = uniprot_handler.getUniProtData(uniprot_identifier).get('ipro',[])
        # Pfam associated to IPRO
        ipro_pf = [name_conv.IPROtoPfam(i) for i in up_ipro]
        # spice up _domains with the additional information
        merged = pfam_domains + [x for x in [item for sublist in ipro_pf for item in sublist] if x not in pfam_domains]
        
        # add fallback domain
        if use_fallback_domains:
            merged.append('FB')
            
        for (i, dom) in enumerate(merged):
            self._domains.append(Domain(dom, self, i))
        
        # set self-interaction initially to false
        self._self_interacting = False
        
        # look for multiplicity of a domain and save mapping to domain instances in self._mult_pfam_map
        pfam_count = defaultdict(int)
        self._mult_pfam_map = defaultdict(list)
        for dom in self._domains:
            pfam_count[dom._pfam_id] += 1
            self._mult_pfam_map[dom._pfam_id].append(dom)
        for pfam_id in pfam_count:
            if pfam_count[pfam_id] == 1:
                del self._mult_pfam_map[pfam_id]
        #print self._mult_pfam_map
    
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
    def __init__(self, protein_ids, ppi, threshold = 0.0, use_fallback_domains = True):
        self._proteins = []
        self._ppi = ppi
        self._up_map = {}
        self._edge_grow_hash = set()
        self._domain_grow_hash = set()
        self._amb_remove = 0
        
        for prot in protein_ids:
            self._proteins.append(Protein(prot, use_fallback_domains))
        
        # determine domain-locations
        domain_map = defaultdict(list)
        for prot in self._proteins:
            for domain in prot._domains:
                domain_map[domain._pfam_id].append(domain)
        
        # set intitial interactions
        for domain_name in domain_map:
            for target_domain_name in domain_query.getDDIs(domain_name):
                
                # no directionality/double occurance, no partner
                if target_domain_name < domain_name or target_domain_name not in domain_map:
                    continue
                
                # add interactions
                for dom1 in domain_map[domain_name]:
                    for dom2 in domain_map[target_domain_name]:
                        #self-interaction
                        if dom1._protein._up_id == dom2._protein._up_id:
                            dom1._protein._self_interacting = True
                            continue
                        # no double edges
                        if target_domain_name == domain_name and dom1._protein._up_id < dom2._protein._up_id:
                            continue
                        # filter for known interaction (at least to a certain degree known or assumed) up to given threshold
                        if dom2._protein._up_id in self._ppi._partners[dom1._protein._up_id] and self._ppi._weights[(dom2._protein._up_id, dom1._protein._up_id)] > threshold:
                            dom1._partners.append(dom2)
                            dom2._partners.append(dom1)
                            
        del domain_map
        # done with connection-building
        
        # remove disconnected proteins from network
        self.cleanup()
        for prot in self._proteins:
            self._up_map[prot._up_id] = prot
        
        # test this ...
        #self.sparsifyFilter()
        
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
                    if P_n >= 0.5:
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
                        if P_n >= 0.5:
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
            elif len(proteins) >= cutoff:
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
    
    def growPairs(self, protein, tfs, cutoff, threshold = 0.75, verbose = False):
        '''
        First builds candidate starting seeds, then seeds.
        DDI-information for starting clusters is partially neglected by allowing connections without domains!
        '''
        jobs = []
        w = self._ppi._weights
        
        # needed to determine if in DDIN
        ref_ddi = self.SafeGetIncidentNodes([protein])
        
        if protein in tfs:
            to_del = []
            for p2 in ref_ddi:
                if p2 in tfs and p2<protein:
                    to_del.append(p2)
            for p2 in to_del:
                del ref_ddi[p2]
            del to_del
            
        tfs = None
        
        # there are DDI with this protein
        if len(ref_ddi) > 0:
            # determine max/above threshold
            # all ddi partners
            ddi_partners = ref_ddi.keys()
            
            # pre_sort all DDII-partners to gather some information
            ddi_partners.sort(key=lambda x:w[(protein, x)], reverse = True)
            
            # if even max below threshold, compute 2 variants anyway
            if w[(protein, ddi_partners[0])] < threshold:
                ddi_partners = ddi_partners[:min(len(ddi_partners), 2)]
            else:
                ddi_partners = [x for x in ddi_partners if w[(x, protein)] > threshold]
                
            with_ddi = {}
            for ddi_partner in ddi_partners:
                with_ddi[ddi_partner] = ref_ddi[ddi_partner]
            
            del ddi_partners
            
            alternatives = [(x, self.filterAlternatives(with_ddi[x])) for x in with_ddi]
                
            with_ddi = []
            for partner, variants in alternatives:
                (w_in, w_out) = self._ppi.computeClusterCohesiveness([protein, partner])
                with_ddi.extend([(set([protein, partner]),set(PSort(d)), w_in, w_out, self._ppi._weights[(protein,partner)]) for d in variants])
            
            jobs = with_ddi
        
        # no domain -> get from PPI only, should only happen without fallback enabled
        else:
            # all ppi partners
            ppi_partners = list(self._ppi._partners[protein])
            
            # pre_sort all PPI-partners to gather some information
            ppi_partners.sort(key=lambda x:w[(protein, x)], reverse = True)
            
            # if even max below threshold, compute 2 variants anyway
            if len(ppi_partners) == 0:
                return set()
            
            if w[(protein, ppi_partners[0])] < threshold:
                ppi_partners = ppi_partners[:min(len(ppi_partners), 2)]
                
            else:
                ppi_partners = [x for x in ppi_partners if w[(x, protein)] > threshold]
                
            jobs = [(set([protein, partner]), set(), self._ppi.computeClusterCohesiveness([protein, partner])[0], self._ppi.computeClusterCohesiveness([protein, partner])[1], self._ppi._weights[(protein,partner)]) for partner in ppi_partners]
        
        results = set()
        no_p = len(jobs)
        print 'Seed from', no_p,'pair(s).'
        for i, job in enumerate(jobs):
            result = self.growManager(cutoff, job)
            results.update(result)
            print 'Finished',str(i+1)+'/'+str(no_p) +'.'
            if verbose:
                print result
            sys.stdout.flush()
            
        return results
    
    def sparsifyFilter(self):
        '''
        If a protein has several identical Pfam-domains and these interact with several identical 
        Pfam-domains of another protein, the subnet of DDIs can be simplified without losing any information
        but preventing combinatorial explosion.
        The edge have to be equally distributed.
        '''
        for i in range(len(self._proteins)-1):
            prot1 = self._proteins[i]
            map1 = prot1._mult_pfam_map
            if len(map1) == 0:
                continue
            for j in range(i+1, len(self._proteins)):
                prot2 = self._proteins[j]
                map2 = prot2._mult_pfam_map
                if prot2._up_id not in self._ppi._partners[prot1._up_id] or len(map2) == 0:
                    continue
                
                # sparsify
                for pfam1 in map1:
                    for pfam2 in map2:
                        # check if interacting
                        domains1 = map1[pfam1]
                        domains2 = map2[pfam2]
                        # since we have a bipartite setting here, this will hold
                        if domains2[0] not in domains1[0]._partners:
                            continue
                        
                        m = len(domains1)
                        n = float(len(domains2))
                        ratio = m/n
                        if (ratio < 1):
                            temp = domains1
                            domains1 = domains2
                            domains2 = temp
                            ratio = 1.0/ratio
                        #print 'sparsify:'+str(domains1)+' to '+str(domains2)
                        #print ratio
                        
                        # must be equally distributed
                        if int(ratio) != ratio:
                            continue
                        ratio = int(ratio)
                        
                        # ratio -> ratio of first/second
                        # clear all
                        for dom1 in domains1:
                            for dom2 in domains2:
                                dom1._partners.remove(dom2)
                                dom2._partners.remove(dom1)
                                
                        # set new in a sparse way
                        for i in range(0, len(domains2)):
                            dom2 = domains2[i]
                            for j in range(i*ratio,(i+1)*ratio):
                                dom1 = domains1[j]
                                
                                dom1._partners.append(dom2)
                                dom2._partners.append(dom1)