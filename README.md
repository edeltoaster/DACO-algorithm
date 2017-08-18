-----------------------------------------------------------------------------------
            DACO: DOMAIN-AWARE COHESIVENESS OPTIMIZATION ALGORITHM
			(prototype by Thorsten Will, 2014; update of May 2016)
-----------------------------------------------------------------------------------

USAGE:

python DACO.py [PPIN] [TFs] [PAIR-THRESHOLD] [DEPTH] [OUT-FILE] , with:

 PPIN : a protein-protein interaction network in the simple input format (SIF)
        with the layout 'ProteinA   ProteinB    weight'. The weight is expected
        to be represented as the probability of an interaction (float between 0 and 1).
        Proteins are expected to be UniProt accession numbers, like 'P12345'.

 TFs : a textfile with a list of transcription factors of the organism line by line.

 PAIR-THRESHOLD : Probability threshold for initial pair seed construction.

 DEPTH : maximal complex size (to keep combinatorial explosion low).

 OUT-FILE : output will be written to this file. Complex candidates are stored line by line,
			members are comma-separated.


Examples of input files for yeast can be found in the 'example_inputs/' folder.

Yeast data of the Dec. 2013 release of UniProt (and Pfam/InterPro at that time) is precached to speed-up the data retrieval.
Without this preloading, the domain-domain interaction network construction can take some time (depending on the size of the network).

To reset the precached yeast data, just delete the following files in the 'data/' folder:
    ipro_mapping.data - relating InterPro and Pfam IDs
    pfam.data - associated Pfam IDs per protein
    up_entry_data.data - Uniprot data per protein
    up_names.data - UniProt and name conversion data

Updates/Fixes:
    - URL of Pfam changed
