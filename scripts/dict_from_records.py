# -*- coding: utf-8 -*-
"""
Created on Thu Jan 22 08:14:44 2015

Reads in the pickled docLists and combines them, then creates a dictionary.
The argument is the number of words to include in the dictionary.  
This will include the top W most frequent words in the corpus

By default this removes stopwords using the NLTK stopword corpus

@author: jmf
"""

import sys
import cPickle
from helper_funcs   import  *
from os.path        import  isfile
from os.path        import  join
from os             import  listdir

cutoff      = sys.argv[1]
inp         = sys.argv[2]
docFiles    = listdir(inp)
 
docList     = {}

"""Assumes one file per doc"""
if len(docFiles) > 1:
    for fi in docFiles:
        print fi
        if isfile(inp+fi):
            with open(inp+fi,'rb') as f:
                d           = f.read()
                print d
                docList[fi] = d
    


else:
    """Assumes \n separated docs in one txt"""    
    docIndex    = 0
    with open(join(inp,docFiles[0]),'rb') as f:           
        for line in f:
            docIndex +=1 
            if len(line) > 1:
                docList[docIndex] = line                          

whole_str   =   ""
for k,doc in docList.iteritems():
    whole_str+=doc.lower()
    
print whole_str    
    
print "before"
file_list_to_lda(whole_str,int(cutoff),stop="yes")
print "after"
    

