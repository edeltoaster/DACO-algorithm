#! /usr/bin/env python
'''
Used to access UniProt-data of all kind
'''

import urllib2, os, cPickle, atexit, name_conv

d_path = '/'.join(os.path.dirname(os.path.abspath(__file__)).split('/')[:-1])+'/data/'

# cached data
__up_data = {}
__up_data_entries_changed = False
__version = 1 # enables precaching and updating over time
if os.path.isfile(d_path+'up_entry_data.data'):
    __up_data = cPickle.load(open(d_path+'up_entry_data.data','rb'))
     
def _setChangedFlag():
    '''
    Helper flag for caching
    '''
    global __up_data_entries_changed
    __up_data_entries_changed = True
    
@atexit.register
def _saveUPData():
    '''
    Saves UniProt-queries so far ...
    '''
    if not __up_data_entries_changed:
        return
    
    f = open(d_path+'up_entry_data.data','wb')
    cPickle.dump(__up_data, f, True)
    
def _readUniProt(up_id):
    '''
    Retrieve and parse UniProt-data for information about some UniProtID
    '''
    try:
        # fields in version 1
        mass = 0
        pfam = []
        ipro = []
        locus = []
        orf = []
        GO_location = []
        GO_function = []
        GO_process = []
        synonyms = []
        name = []
        id_field = ''
        AAs = 0
        
        request = urllib2.urlopen('http://www.uniprot.org/uniprot/'+up_id+'.txt')
        data = request.read()
        request.close()
        
        for line in data.split('\n'):
            temp = line.split()
            
            if len(temp) == 0:
                continue
                
            if temp[0][0] == 'R':
                continue
            elif temp[0] == 'SQ':
                mass = float(temp[4])
                break
            elif temp[0] == 'ID':
                id_field = temp[1]
                AAs = int(temp[3])
                continue
            elif temp[0] == 'DR':
                if temp[1] == 'InterPro;':
                    appearance = 1
                    if len(temp) == 5:
                        appearance = int(temp[4][:-1])
                    ipro.extend(appearance * [temp[2][:-1]])
                    
                elif temp[1] == 'Pfam;':
                    appearance = 1
                    if len(temp) == 5:
                        appearance = int(temp[4][:-1])
                    pfam.extend(appearance * [temp[2][:-1]])
                
                # GO annotation
                elif temp[1] == 'GO;':
                    if temp[3][0] == 'C':
                        GO_location.append(temp[2][:-1])
                    elif temp[3][0] == 'F':
                        GO_function.append(temp[2][:-1])
                    elif temp[3][0] == 'P':
                        GO_process.append(temp[2][:-1])
            elif temp[0] == 'GN':
                for entry in ' '.join(temp[1:]).split(';'):
                    temp2 = entry.strip().split('=')
                    if temp2[0] == 'Name':
                        name.append(temp2[1])
                    elif temp2[0] == 'OrderedLocusNames':
                        locus.append(temp2[1])
                    elif temp2[0] == 'ORFNames':
                        orf.append(temp2[1])
                    elif temp2[0] == 'Synonyms':
                        synonyms = temp2[1].split(', ')
                        
        # save data to dictionary
        data = {}
        data['version'] = __version
        data['mass'] = mass
        data['pfam'] = pfam
        data['ipro'] = ipro
        data['locus'] = locus
        data['orf'] = orf
        data['GO_location'] = GO_location
        data['GO_process'] = GO_process
        data['GO_function'] = GO_function
        data['synonyms'] = [i.upper() for i in synonyms]
        data['id'] = id_field
        data['AAs'] = AAs
        data['name'] = [iname.upper() for iname in name]
        
        # only non-ambigous ones ...
        for syn in data['orf']+data['locus']:
            name_conv.__name_data[syn] = up_id
        name_conv.__name_data[id_field] = up_id
        name_conv.__name_data[up_id] = id_field
        
        __up_data[up_id] = data
        _setChangedFlag()
        #enforce update
        name_conv.__name_data_entries = -1
        
        return data 
    
    except urllib2.HTTPError, err:
        # save data to dictionary
        data = {}
        data['version'] = __version
        data['mass'] = mass
        data['pfam'] = pfam
        data['ipro'] = ipro
        data['locus'] = locus
        data['orf'] = orf
        data['GO_location'] = GO_location
        data['GO_process'] = GO_process
        data['GO_function'] = GO_function
        data['synonyms'] = [i.upper() for i in synonyms]
        data['id'] = id_field
        data['AAs'] = AAs
        data['name'] = [iname.upper() for iname in name]
        name_conv.__name_data[id_field] = up_id
        name_conv.__name_data[up_id] = id_field
        #enforce update
        name_conv.__name_data_entries = -1
        if err.code == 404:
            __up_data[up_id] = data
            _setChangedFlag()
            return data
        else:
            __up_data[up_id] = data
            _setChangedFlag()
            return data
        
    except Exception:
        data = {}
        data['version'] = __version
        data['mass'] = mass
        data['pfam'] = pfam
        data['ipro'] = ipro
        data['locus'] = locus
        data['orf'] = orf
        data['GO_location'] = GO_location
        data['GO_process'] = GO_process
        data['GO_function'] = GO_function
        data['synonyms'] = [i.upper() for i in synonyms]
        data['id'] = id_field
        data['AAs'] = AAs
        data['name'] = [iname.upper() for iname in name]
        # return empty data
        __up_data[up_id] = data
        _setChangedFlag()
        name_conv.__name_data[id_field] = up_id
        name_conv.__name_data[up_id] = id_field
        #enforce update
        name_conv.__name_data_entries = -1
        
        return data

def getUniProtData(up_id):
    '''
    Get UniProt-data for information about some UniProtID, output as map.
    Currently: version, id, name, AAs, mass, pfam, ipro, locus, orf, GO_location, GO_process, GO_function, synonyms
    '''
    if up_id in __up_data:
        # if older internal __version: update
        if __up_data[up_id]['version'] == __version:
            return __up_data[up_id]
           
    return _readUniProt(up_id)

def getTrivialNames(up_id):
    if not up_id in __up_data:
        _readUniProt(up_id)
    
    names = __up_data[up_id]['name']+__up_data[up_id]['synonyms']
    if names == []:
        names = [name_conv.NameFromUniprot(up_id)]
    
    return [x.split('_')[0] for x in names]

def parseUniProtFile(f , verbose = False):
    '''
    Given a UniProt-file f, parse the file and return a list of UP-IDs
    '''
    slashcount = 0
    up_ids = []
    n_total = 0
    n_struct = 0
    for line in open(f):
        #print slashcount
        if slashcount < 2:
            if line.startswith('_'):
                slashcount += 1
            continue
        if line.startswith('-'):
            break
        line = line.strip()
        up_ids.append(line[-42:-36])
        n_total += 1
        if line.split()[-2] == '(3)':
            n_struct += 1
    if verbose:
        print 'total:', n_total
        print 'with structure:', n_struct
        print 'perc:', n_struct / float(n_total)
    
    return up_ids