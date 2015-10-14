# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 14:54:40 2015

input   -i should be a folder containing your documents
output  -o sould be a filename for the pickled docList

docList will be a dictionary if "document_title" : "stemmed_document_text"

@author: frickjm
"""

import sys
import getopt
import helper_funcs
import cPickle
from os import listdir
from os.path import join
from os.path import isfile


def main(argv):
    inp, ofile  = handleArgs(argv)
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

    #out     = helper_funcs.stem_doc_dict(docList)
    out     = docList
    with open(ofile,'wb') as f:
        cp  = cPickle.Pickler(f)
        cp.dump(out)
    
    
def handleArgs(argv):
    try:
       opts, args = getopt.getopt(argv,"hi:o:",["inp=","outp="])
    except getopt.GetoptError:
       print 'test.py -i <input folder> -o <output pickle>'
       sys.exit(2)
    for opt, arg in opts:
       if opt == '-h':
          print 'pickleDocs.py -i <inputfile> -o <outputfile>'
          print 'example: pickleDocs.py -i ../data/in -o ../data/pickledDocs/wikiDocs.pickle'
          sys.exit()
       elif opt in ("-i", "--ifile"):
          inp = arg
       elif opt in ("-o", "--ofile"):
          ofile = arg
    print 'Input folder is ', inp
    print 'Output file is ', ofile 
    return inp, ofile


if __name__ == "__main__":
    main(sys.argv[1:])