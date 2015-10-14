# -*- coding: utf-8 -*-
"""
@author: James Frick

NIH Center for Advancing Translational Sciences
jamesmfrick@gmail.com

Evaluates a model.  Can be used to create topic heatmaps, histograms, or score
a model based on Disease Ontology hierarchy ROC

Argument 1:
    if this is an integer, script will score the PCA/word frequency version,
        using this arg as the first principal component
    else this should be the folder that contains your data
    
Argument 2:
    if doing PCA/word freq, this is the last PC used
    else it should be either:
        1) "runforall" which will perform analysis on all OMIM-DO docs
        2) the name of a disease in which you're interested.
            if this is the case, it will find OMIM entries containing that name
            and ask you which you'd like to search on.
            
The heatmap portion has been entirely commented out in order to avoid reliance
on Matplotlib.  Thus models can be evaluated on machines with fewer libraries.
"""

import sys, cPickle
import helper_funcs
from os import listdir
import matplotlib.pyplot as plt
import numpy as np
import csv
from textwrap import wrap

def main():
 
    #f2  =   open("../data/out/doc_topics.txt",'w')
    
    ###########################################################################   
    """READ INPUT PARAMETERS
    #input parameters are :
        #vocab is the dictionary
        #lambda is a matrix of the words for each Topic
        #gamma is a matrix of the topics for each Document
        #search is a term you would like to match in a document title
            #if search contains spaces, it will search on two terms together
            #and by default should sum the two topic vectors"""
    if len(sys.argv) > 1:
        folder     = sys.argv[1]
        allsearch  = []
        i1          = folder.find("models/")+7
        i2          = folder.find("kwords")
        dsize       = folder[i1:i2]

        vn          = "../data/dictionary.txt"
        vocab       = str.split(file(vn).read())
        lambdap    = folder + "lambda.pickle"
        gammap      = folder + "gamma.pickle"
        for i in range(2,len(sys.argv)):
            allsearch.append(sys.argv[i])
                    
    else:
        quit

    with open(lambdap,'rb') as tl:
        testlambda  = cPickle.load(tl)
        
    countDict = {}        
    with open("../data/dictCounts.txt",'rb') as f:
        for line in f:
            sp  = line.split("\t")
            if len(sp) > 1:
                countDict[sp[0]] = int(sp[1].strip())

    stemmedMap  = {}
    with open("../data/stemmed_mapping.txt",'rb') as f:
        for line in f:
            sp = line.split("\t")
            if len(sp) > 1:
                stemmedMap[sp[0]] = sp[1].strip()

        ###########################################################################
    """Write out the word distributions for each topic, sorted by the probability
        that a given word 'belongs' to that topic"""
    copyLambda  = testlambda.copy()
 
    stemmedMap  = {}
    with open("../data/stemmed_mapping.txt",'rb') as f:
        for line in f:
            sp = line.split("\t")
            if len(sp) > 1:
                print sp
                stemmedMap[sp[0]] = sp[1].strip()
                try:
                    vocab[vocab.index(sp[0])] = sp[1].strip()
                except ValueError as e:
                    print e
                    
    countDict   = {}                    
    with open("../data/dictCounts.txt",'rb') as f:
        for line in f:
            sp = line.split("\t")
            if len(sp) > 1:
                print sp
                countDict[stemmedMap[sp[0]]] = int(sp[1].strip())
                
    
#    with open(folder+"topics2.txt",'w') as f:
##        for k1 in range(0,len(testlambda[0])):
##            lambdak1 = testlambda[:,k1]
##            lambdak1 = lambdak1/sum(lambdak1)
##            testlambda[:,k1] = lambdak1
#        for k in range(0, len(testlambda)):
#            lambdak = testlambda[k, :]
#            for l in range(0,len(lambdak)):
#                print countDict[vocab[l]]
#                lambdak[l] = lambdak[l]*(1 + np.log(countDict[vocab[l]]/10 +0.0001))
#            #lambdak = lambdak / sum(lambdak)
#            temp = zip(lambdak, range(0, len(lambdak)))
#            temp = sorted(temp, key = lambda x: x[0], reverse=True)
#            f.write("topic "+str(k)+"\n")
#            
#            #Write out the top 100 words per topic
#            for i in range(0, min(50,len(vocab))):
##                print str(vocab[temp[i][1]])
##                if stemmedMap.has_key(vocab[temp[i][1]]):
##                    print str(stemmedMap[vocab[temp[i][1]]])
##                print str(countDict[vocab[temp[i][1]]])
##                print temp[i][0], float(temp[i][0])*countDict[vocab[temp[i][1]]]
##                stop=raw_input("")
#                f.write(str(vocab[temp[i][1]]) +"\t\t\t" + str(temp[i][0])+"\n")
#            f.write("\n\n")    
#    
#    
#    with open(folder+"topics3.txt",'w') as f:
#        for k1 in range(0,len(copyLambda)):
#            lambdak1 = copyLambda[k1,:]
#            lambdak1 = lambdak1/sum(lambdak1)
#            copyLambda[k1,:] = lambdak1
#        for k1 in range(0,len(copyLambda[0])):
#            lambdak1 = copyLambda[:,k1]
#            lambdak1 = lambdak1/sum(lambdak1)
#            copyLambda[:,k1] = lambdak1
#        for k in range(0, len(copyLambda)):
#            lambdak = copyLambda[k, :]
#            #lambdak = lambdak / sum(lambdak)
#            temp = zip(lambdak, range(0, len(lambdak)))
#            temp = sorted(temp, key = lambda x: x[0], reverse=True)
#            f.write("topic "+str(k)+"\n")
#            
#            #Write out the top 100 words per topic
#            for i in range(0, min(50,len(vocab))):
#                f.write(str(vocab[temp[i][1]]) +"\t\t\t" + str(temp[i][0])+"\n")
#            f.write("\n\n")

    THRESHOLD   = 10   
    THRESHOLD2  = .3        
    with open(folder+"topics.txt",'w') as f:
        for k1 in range(0,len(testlambda[0])):
            lambdak1 = testlambda[:,k1]
            lambdak1 = lambdak1/sum(lambdak1)
            testlambda[:,k1] = lambdak1
        testsum=0
        for k in range(0, len(testlambda)):
            lambdak = testlambda[k, :]
            testsum+=lambdak[0]
            print testsum
            #lambdak = lambdak / sum(lambdak)
            temp = zip(lambdak, range(0, len(lambdak)))
            temp = sorted(temp, key = lambda x: x[0], reverse=True)
            f.write("topic "+str(k)+"\n")
            
            #Write out the top 100 words per topic
            for i in range(0, min(500,len(vocab))):
                if countDict[vocab[temp[i][1]]] > THRESHOLD:
                    if temp[i][0] > THRESHOLD2:
                        
                        f.write(str(vocab[temp[i][1]]) +"\t\t\t" + str(temp[i][0])+"\n")
            f.write("\n\n")
#            
    #    for l in range(0, len(testgamma)):
    #        lambdal = list(testgamma[l, :])
    #        #print lambdal
    #        f2.write("document "+str(l)+":\n"+str(lambdal)+"\n\n")


if __name__ == '__main__':
    main()
