# -*- coding: utf-8 -*-
"""
Created on Thu Jan 22 15:50:57 2015

@author: James Frick

NIH Center for Advancing Translational Sciences
jamesmfrick@gmail.com

This is the script that trains a new LDA model on the documents set.

By default, it fits a model with 10 topics and uses the dictionary_7k.txt as 
its dictionary.  New dictionaries can be trained with dict_from_records.py.

To change topic number or dictionary, add arguments when calling the script:

python callLDA.py 40 ../data/dictionary_15k.txt

The above call fits a model with 40 topics, using the dictionary specified.


To also specify the hyperparameters alpha and beta (as, say, 0.1 and 1.0):

python callLDA.py 40 ../data/dictionary_15k.txt 0.1 1.0


"""

import cPickle
import helper_funcs
import onlineldavb
import sys
from os import listdir
from os.path import isdir
from os import mkdir
import operator
import numpy as np
import getopt





def handleArgs(argv):
    try:
       opts, args = getopt.getopt(argv,"k:d:a:b:s",["topics=","folder=","alpha=","beta=","save"])
    except getopt.GetoptError:
       print 'test.py -k <# of Topics> -d <dictionary path>'
       sys.exit(2)
       
    alpha   = None
    beta    = None
    save    = None   
    for opt, arg in opts:
        print opt, arg

        if opt in ("-k", "--topics"):
            k = arg
        elif opt in ("-d", "--folder"):
            folder = arg
        elif opt in ("-a", "--alpha"):
            alpha = arg
        elif opt in ("-b", "--beta"):
            beta = arg
        elif opt in ("-s", "--save"):
            save = arg
    print 'Number of topics entered: ', k
    print 'Dictionary is located at: ', folder
    print alpha, beta, save
    return k, folder, alpha, beta, save




def main(argv):
    
    doc_list    =   []
    
    argList     = handleArgs(argv)
    #list the docs in pickledDocs folder
    p   =   "../data/pickledDocs/"
    l   =   listdir(p)
    fileList    =   [p+f for f in l]
    
    #for each pickled doclist, append all docs to master doclist
    for fi in fileList:
        with open(fi,'rb') as d:
            docs    =   cPickle.load(d)
        for k,x in docs.iteritems(): 
            doc_list.append(x)
        print len(doc_list)
        
        
    
    #D is total number of docs to show to the model, K is number of topics
    goal_its    =   80                #number of iterations to run LDA 
    corp_size   =   len(doc_list)       #number of documents in the corpus
    D           =   corp_size*goal_its  #number of documents expected to see
    K           =   10                  #default topic value, if none given in parameters
    saveModel   =   False               #whether to save LDA model itself
    desc        =   ""                  #for performing non-standard runs
    version     =   ""                  #for having multiple models with same parameters
    hyper_param =   ""                  #for testing hyperparameters
    
    #define the vocabulary file we will be using
    vocab       =   helper_funcs.read_dict("../data/dictionary.txt") #default dict 
    
    
    #initialize an instance of the OnlineLDA algorithm
    #parameters - dictionary, num topics, learning rate, beta, tau, kappa
    #if the path to an OnlineLDA pickle is passed, it re-opens that pickle

    K           =   int(argList[0])
    vocab       =   vocab = str.split(file(argList[1]).read())
    if not (argList[2] is None):
        alpha   = float(argList[2])
    else:
        alpha   =   0.1
    if not (argList[3] is None):
        beta    = float(argList[3])
    else:
        beta     =   1.

    saveModel   =   False
    lda         =   onlineldavb.OnlineLDA(vocab,K,D,alpha,beta,1024,0.)
    print "created LDA with parameters:\nnumwords: "+str(len(vocab))+"\n#topics: "+str(K)+"\nalpha: "+str(alpha)+"\nbeta: "+str(beta)
           
    paramTitle  =   hyper_param+str(len(vocab)/1000)+"kwords_"+str(K)+"topics"
    
    folder  = "../data/out/"+paramTitle
    if not isdir(folder):
        mkdir(folder)
    
    W           =   len(vocab)
    
    print "dictionary size: " + str(W)
    print paramTitle
    
    
    
    print folder
    #if desc.find("label") > -1:
    #    with open("../data/out/past_models/"+paramTitle+"/dictionary.txt",'wb') as f:
    #        voc = sorted(vocab.items(),key=operator.itemgetter(1))
    #        for x in voc:
    #            f.write(x[0]+"\n")
    #perform LDA on the document list for goal_its iterations, updating lambda
    for i in range(lda._updatect,goal_its):
        (gamma, bound)      = lda.update_lambda(doc_list)
        
        (wordids, wordcts)  = onlineldavb.parse_doc_list(doc_list,lda._vocab)
        perwordbound        = bound * len(doc_list) / (D*sum(map(sum,wordcts)))
        print i, np.exp(-perwordbound)
        
        #pickle the model and its output occasionally
        if (i+1) == goal_its:
            if not isdir(folder):
                mkdir(folder)
            with open(folder+"/gamma.pickle",'wb') as f:
                cp2 = cPickle.Pickler(f)
                cp2.dump(gamma)
            with open(folder+"/lambda.pickle",'wb') as f:
                cp  = cPickle.Pickler(f)
                cp.dump(lda._lambda)
            np.savetxt(folder+'/lambda.dat', lda._lambda)
            
            
            if saveModel:
                
                with open(folder+"/LDA.pickle",'wb') as f:
                    cp3 = cPickle.Pickler(f)
                    cp3.dump(lda)
            

if __name__ == "__main__":
    main(sys.argv[1:])

