#! /usr/bin/env python
'''
Scripts for protein/domain name conversion, supported:
gene name <-> UniprotID, InterPro -> Pfam
'''

import urllib, urllib2, cPickle, os, atexit, re

d_path = '/'.join(os.path.dirname(os.path.abspath(__file__)).split('/')[:-1])+'/data/'

# cached data
__name_data = {}
__name_data_entries = 0
if os.path.isfile(d_path+'up_names.data'):
    __name_data = cPickle.load(open(d_path+'up_names.data','rb'))
    __name_data_entries = len(__name_data)
    
__ipro_data = {}
__ipro_data_entries = 0
if os.path.isfile(d_path+'ipro_mapping.data'):
    __ipro_data = cPickle.load(open(d_path+'ipro_mapping.data','rb'))
    __ipro_data_entries = len(__ipro_data)
    
@atexit.register
def _saveUPNames():
    '''
    Saves UniProt-queries so far ...
    '''
    if len(__name_data) == __name_data_entries:
        return
    f = open(d_path+'up_names.data','wb')
    cPickle.dump(__name_data, f, True)
    
@atexit.register
def _saveInterProMap():
    '''
    Saves InterPro-mapping-queries so far ...
    '''
    if len(__ipro_data) == __ipro_data_entries:
        return
    f = open(d_path+'ipro_mapping.data','wb')
    cPickle.dump(__ipro_data, f, True)
                               
def NameToUniprot(query):
    '''
    Converts gene name to UniProt-ID if possible, returns '' otherwise
    '''
    if query in __name_data:
        return __name_data[query]
    
    params = {'from':'ACC+ID', 'to':'ACC', 'format':'tab','query':query}
    return _getFromUniprot(params)[query]

def NameFromUniprot(query):
    '''
    Converts UniProt-ID to gene name if possible, returns '' otherwise
    '''
    if query in __name_data:
        return __name_data[query]
    
    params = {'from':'ACC', 'to':'ID', 'format':'tab','query':query}
    return _getFromUniprot(params)[query]

def NamesToUniprot(mquery):
    '''
    Maps gene names to UniProt-IDs if possible, returns '' otherwise.
    '''
    params = {'from':'ACC+ID', 'to':'ACC', 'format':'tab','query':' '.join(mquery)}
    return _getFromUniprot(params)

def NamesFromUniprot(mquery):
    '''
    Maps UniProt-IDs to gene names if possible, returns '' otherwise.
    '''
    params = {'from':'ACC', 'to':'ID', 'format':'tab','query':' '.join(mquery)}
    return _getFromUniprot(params)
    
def EntrezGenesToUniprot(mquery):
    '''
    Maps entrezgene ids to UniProt-IDs if possible, returns '' otherwise.
    '''
    params = {'from':'P_ENTREZGENEID', 'to':'ACC', 'format':'tab','query':' '.join(mquery)}
    return _getFromUniprot(params)
    
def _getFromUniprot(params):
    '''
    Sends a mapping-query to UniProt
    '''
    
    results = {}
    
    # check if cached
    new_query = []
    query = params['query'].split(' ');
    for q in query:
        if q.startswith('"') and q.endswith('"'):
            q = q[1:-1]
        if q in __name_data:
            results.update({q:__name_data[q]})
        else:
            new_query.append(q)
        
    # retrieve rest if needed
    if len(new_query) > 0:
        params['query'] = ' '.join(new_query)
        url = 'http://www.uniprot.org/mapping/'
        data = urllib.urlencode(params)
        request = urllib2.Request(url, data)
        request.add_header('User-Agent', 'Python contact')
        try:
            response = urllib2.urlopen(request)
            page = response.read().strip()
            
            for line in page.split('\n'):
                if line.startswith('From'):
                    continue
                temp = line.strip().split('\t')
                results.update({temp[0]:temp[1]})
                __name_data.update({temp[0]:temp[1]})
        except:
            pass
        
        # check for non-retrievable subqueries -> update __name_data accordingly
        if len(query) > len(results):
            for element in query:
                if element not in results:
                    results.update({element:''})
                    __name_data.update({element:''})
    
    return results

def getNameListFromUniProtIDs(upids):
    nm = NamesFromUniprot(upids)
    return [nm[x] for x in upids if x in nm]

def IPROtoPfam(ipro_id):
    '''
    Queries mapping InterPro->Pfam, sometimes several Pfam domains are associated
    '''
    if ipro_id in __ipro_data:
        return __ipro_data[ipro_id]
    
    params = {'db':'interpro', 'id':ipro_id, 'style':'raw'}
    url = 'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch'
    data = urllib.urlencode(params)
    request = urllib2.Request(url, data)
    response = urllib2.urlopen(request)
    page = response.read().strip()
    tables = re.findall(r'<TR>(.*?)</TR>', page, re.M|re.I|re.S)
    
    pfam_ids = []
    try:
        temp = tables[3].strip().split('valign="top">')[1].split('<br>')
        
        for signature in temp:
            if signature.startswith('PFAM:'):
                pf = signature.split('</a>')[0].split('>')[-1]
                pfam_ids.append(pf)
    except Exception:
        pass
    
    # update
    __ipro_data[ipro_id] = pfam_ids
    
    return pfam_ids