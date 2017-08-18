#! /usr/bin/env python
import urllib2, sys, cPickle, os, atexit
import xml.etree.ElementTree as ET

'''
Query PFAM database and obtain DDIs
'''
d_path = '/'.join(os.path.dirname(os.path.abspath(__file__)).split('/')[:-1])+'/data/'

# cached Pfam-data (if cached)
__pfam = {}
__pfam_entries_changed = False
if os.path.isfile(d_path+'pfam.data'):
    __pfam = cPickle.load(open(d_path+'pfam.data','rb'))

# cached DDI data (pre-build from IDDI and DOMINE databases)
__ddis = cPickle.load(open(d_path+'ddi.data','rb'))

def _setPfamChangedFlag():
    '''
    Helper flag for caching
    '''
    global __pfam_entries_changed
    __pfam_entries_changed = True

@atexit.register
def _savePfam():
    '''
    Saves Pfam-queries so far ...
    '''
    if not __pfam_entries_changed:
        return

    f = open(d_path+'pfam.data','wb')
    cPickle.dump(__pfam, f, True)

def getPfamData(up_id, use_mirror = True):
    '''
    Get Pfam-entry for Uniprot ID up_id,
    output: list of Pfam-A matches with [Pfam-id, evalue]
    '''

    # get cached
    if up_id in __pfam:
        return __pfam[up_id]

    try:
        # primary source: http://pfam.sanger.ac.uk/protein/
        # best mirror: http://pfam.janelia.org/protein/
        if use_mirror:
            request = urllib2.urlopen('http://pfam.xfam.org/protein/'+up_id+'?output=xml')
        else:
            request = urllib2.urlopen('http://pfam.xfam.org/protein/'+up_id+'?output=xml')
        entry = list(ET.parse(request).getroot())[0]
    except IndexError:
        return []
    except Exception, e:
        if use_mirror == False:
            return getPfamData(up_id, True)

        print e
        print "Problems while retrieving from Pfam"
        sys.exit(1)

    # get matches element in tree
    matches = None
    for e in entry:
        if e.tag == '{http://pfam.sanger.ac.uk/}matches':
            matches = e
            break

    # probably there are no matches
    if matches == None:
        matches = []

    data = []
    for match in matches:
        if match.attrib['type'] == 'Pfam-A':
            match_data = [match.attrib['accession']]
            # old code: location = match.getchildren()[0]
            for location in list(match):
                match_data.append(float(location.attrib['evalue']))
                data.append(match_data)

    # add to cache
    __pfam.update({up_id:data})
    _setPfamChangedFlag()

    return data

def getDDIs(pfam_id):
    '''
    Returns a set of DDI-partners of a given Pfam-id
    '''
    return __ddis[pfam_id]
