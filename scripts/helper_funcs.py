# -*- coding: utf-8 -*-
"""
Created on Thu Jan 22 16:06:41 2015

@author: James Frick

NIH Center for Advancing Translational Sciences
jamesmfrick@gmail.com

Most of the supporting functions for dictionary creation, preprocessing, and scoring
"""
import numpy as np
import cPickle
import operator
from nltk.corpus import stopwords
from nltk.stem  import porter
import re
import enchant
import operator
import sys
import scipy.stats as sts
from scipy.spatial.distance import euclidean

#reads in a dictionary file that is output by createdict.py, adds ID string
def read_dict(filename):
    f       =   open(filename).readlines()
    count   = 0
    diction = {}
    for line in f:
        line = line.replace("\n",'')
        diction[line]   = count
        count+=1
    return diction
    
#read records from the the omim text file
def read_omim_recs(filename):
    mystr       = open(filename).read()
    doc_list    = []
    i2          = 0
    while i2 > -1:
        i1 = mystr.find("*RECORD",i2)
        try:    
            i2 = mystr.find("*RECORD",i1+1)
            rec= mystr[i1:i2]
            doc_list.append(rec)
        except ValueError as err:
            print err
            i2 = -1
    return doc_list
    
def read_in_wiki_pickle(pickleString):
    with open(pickleString,'rb') as titleP:
        wikiDocs    =   cPickle.load(titleP)
    return wikiDocs



def read_topic_dist(gammaPickle,titlePickle):
    docNames    =   []
    for t in titlePickle:
        with open(t,'rb') as titleP:
            dn    =   cPickle.load(titleP)
        for d in dn:
            docNames.append(d)
            
    with open(gammaPickle,'rb') as gammaP:
        docTopics   =   cPickle.load(gammaP)

    return docNames, docTopics

def getAncestorMap(omim,d2d):
    pass

###############################################################################
#SIMILARITY FUNCTIONS    
###############################################################################

    
def find_index(nameString,docNames):
    nameString = nameString.lower()
    print nameString
    count=0
    for i in docNames:
        i = i.lower()
        if i.find(nameString) > -1:
            print i, count
            yn  = str(raw_input("Is this the record you'd like to match? (y/n)\n").lower())
            if yn == "y":
                return count
        count+=1        
    return -1
    
def findMultIndex(allsearch,docNames):
    indices =   []
    for i in range(0,len(allsearch)):
        indices.append(find_index(allsearch[i],docNames))
    return indices
            
def get_sim_docs(docIndex,docTopics,docNames,KL=None,fc=False):
    #print docIndex
    if docIndex == -1:
        print "No Titles Contain Query"
        return []
    else:
        vec1    =   np.zeros(docTopics[0].shape)
        for i in range(0,len(docIndex)):  
            vec1    =   np.add(vec1,docTopics[docIndex[i]])
        #vec1    = np.divide(vec1,np.sum(vec1))
        #vec1    = np.log10(vec1)
    scores  =   []

    l   =   len(docTopics)
    if type(docIndex) is int:
            docIndex=[docIndex]
    for i in range(0,l):
        if i not in docIndex:
            vec2    =   np.array(docTopics[i],dtype=np.float)
            #vec2    =   np.divide(vec2,np.sum(vec2))
            #vec2    =   np.log10(vec2)
            if (not KL) | fc:            
                inner   =   np.dot(vec1,vec2)
                norm1   =   np.sqrt( np.dot(vec1,vec1) )
                norm2   =   np.sqrt( np.dot(vec2,vec2) )
                if fc:
                    scores.append(-euclidean(vec1,vec2))
#                    print vec1
#                    print vec2
#                    print -euclidean(vec1,vec2)
#                    staop=raw_input("wait")
                else:
                    scores.append(np.divide(inner,norm1*norm2))
            else:
                scores.append(-sts.entropy(vec1,vec2))
        else:
            scores.append(99999999)
    names={}
    count=0
    removenum       = re.compile("_#")
    for x in scores:
        d   = docNames[count]
        y   =  re.sub(removenum,'',d)
        source = d[:d.find("_")]
        names[y]    =   [x, list(docTopics[count]), source]
        #print names[docNames[count]]
#        tops[str(x)]              =   docTopics[count]
        count+=1
    nl= sorted(names.items(),key=operator.itemgetter(1),reverse=True)
    return nl

    
    
def threshold(docTopics,threshold):
    y = list(np.where(docTopics>threshold,1.0,0.0000001))
    #print y
    return y
    

def normalizeTopics(docTopics):
    m   =   np.mean(docTopics,axis=0)
    for row in docTopics:
        row =   np.divide(row,m)
    return docTopics
    

#for each result return a tuple of (searchParent,resultParent) or (-1,-1) for no match    
def scoreResults(search,results,ancestryMap):
    if not search is list:
        search=[search]
    ancestry = []
    
    for x in search:
        x=x[:10]            
        if ancestryMap.has_key(x):
            ancestry = ancestryMap[x]
    

    resultMatches=[]
    count_in_results    = 0
    location_in_results = [] 
    for x in results:
        x=x[:10] 
        result_ancestry=[]
        if ancestryMap.has_key(x):
            result_ancestry = ancestryMap[x]
        match=False

        for ra in result_ancestry:
            for a in ancestry:
                if ra == a:
                    match=True
                    
        if not result_ancestry == []:
            if match:
                resultMatches.append(1.)
            else:
                resultMatches.append(0.)
            location_in_results.append(count_in_results)
        count_in_results += 1
                
    return resultMatches, location_in_results
                    
def compareNewDoc(docGamma,docTopics,docNames):
    num = len(docTopics[0])
    vec1 = np.reshape(np.array(docGamma),(num))
    l= len(docTopics)
    scores=[]
    for i in range(0,l):
        vec2    =   np.array(docTopics[i],dtype=np.float)
        scores.append(-sts.entropy(vec1,vec2))
        
    count=0
    names={}
    removenum       = re.compile("_#")
    for x in scores:
        d   = docNames[count]
        y   =  re.sub(removenum,'',d)
        source = d[:d.find("_")]
        names[y]    =   [x, list(docTopics[count]), source]
        count+=1
    nl= sorted(names.items(),key=operator.itemgetter(1),reverse=True)#[0:200]
    return nl

def topDocsByTopic(docTopics,docNames):
    l = len(docTopics[0])
    n = len(docTopics)
    tops = np.array(docTopics)
    inds = [x for x in range(0,n)]
    inds = np.array(inds)
    inds = np.reshape(inds,(inds.shape[0],1))
    print inds.shape, tops.shape
    topics = np.append(tops,inds,axis=1)

    print topics
    a=raw_input("'")
    with open("model/topTopicDocs.txt","wb") as f:
        f.write("Top 10 Documents per Topic:\n\n")
        for i in range(0,l):
            f.write("\nTopic "+str(i)+":\n")
            t = topics[:,i]
            out = np.argsort(t)[-10:][::-1]
            print out
            for j in range(0,10):
                f.write(docNames[out[j]]+"\n")
                print topics[out[j]]
            a=raw_input(str(i))
        

        
    
###############################################################################
#PREPROCESSING FUNCTIONS    
###############################################################################
def get_dictionary(docs, dictFilters=None):
    #dictFilters = None
    dicts = []
    allwords = set()
    for doc in docs:
        allwords = allwords | set(doc.keys())
    if dictFilters:
        for filter in dictFilters:
            allwords = filter(allwords, dicts)
    dict_list = list()
    index_map = {}
    i = 0
    for word in allwords:
        dict_list.append(word)
        index_map[word] = i
        i += 1

    return dict_list, index_map

def write_dictionary(dict_list,totalwords, path="../data/dictionary.txt"):
    count=0
    tag = "_" + str(int(totalwords/1000))+"k.txt"
    path.replace(".txt",tag)
    with open(path, "w") as f:
        for word in dict_list:
            if count < totalwords:
                f.write(word + "\n")
                count+=1
                
def write_stem_mapping(stemmed):
    with open("../data/stemmed_mapping.txt",'wb') as f:
        for k,v in stemmed.iteritems():
            f.write(k+"\t"+v+"\n")
            
    


def plaintext_to_wordcounts(words,cutoff):

    counts  = {}
    #stemmer = nltk.stem.porter.PorterStemmer()
    ret_dict= {}
    for word in words:
        #word   = stemmer.stem(word)
        try:
            counts[word] += 1
        except KeyError:
            counts[word] = 1
    ret     = sorted(counts.items(),key=operator.itemgetter(1),reverse=True)
    print ret
    cutoff  = min(cutoff,len(ret)-1)
    with open("../data/dictCounts.txt",'wb') as f:
        for i in range(0,cutoff):
            ret_dict[ret[i][0]] = ret[i][1]
            f.write(str(ret[i][0])+"\t"+str(ret[i][1])+"\n")
    return ret_dict

def file_list_to_lda(whole_str, cutoff, stop=None):
    ignore  =   []
                
    docs = []
    if stop:
        stops   =   stopwords.words("english")
    engDict = enchant.Dict("en_US")
    stemmer = porter.PorterStemmer()
    no_nums = re.compile("[0-9]+")
    noa1    = re.compile("[a-zA-Z][0-9]")
    words = re.findall("[a-zA-Z0-9]+", whole_str)
    print len(words)
    words = [word for word in words if (word not in stops) & (word not in ignore) &(len(word) > 2)]
    print len(words)
    words   =[word for word in words if not (re.match(no_nums,word) or re.match(noa1,word))]
    print len(words)
    stemmed_dict={}
    for i in range(0,len(words)):
        w   =   words[i]+" ".strip()
        #if engDict.check(w):
        words[i]= stemmer.stem(w)
        try:            
            st      = stemmed_dict[words[i]]            
            if st.find(w) < 0:
                stemmed_dict[words[i]]=st+","+w
        except KeyError:
            stemmed_dict[words[i]]=w
    doc = plaintext_to_wordcounts(words,cutoff)
    docs.append(doc)
    # 'docs' is a list of lists of words that have been pre-processed
    dict_list, index_map = get_dictionary(docs)
    write_dictionary(dict_list,cutoff)
    write_stem_mapping(stemmed_dict)
    return 
    
def stem_docs(docs):
    c=re.compile("[\[\]\(\)\]!,.=\-;']")
    stemmer =   porter.PorterStemmer()
    engDict =   enchant.Dict("en_US")
    alphan = re.compile("^[a-zA-Z0-9]+")
    count=0
    for x in docs:
        x = re.sub(c,' ',x)
        x2  = ""
        ws  =   re.split("\s+",x)
        for w in ws:
            if not w == '':
                try:
                    if not re.match(alphan,w):
                        pass
                    else:
                        if engDict.check(w):
                            w2  =   stemmer.stem(w)
                            x2  += " "+w2
                        else:
                            x2  += " "+w
                except UnicodeDecodeError:
                    pass
        docs[count] = x2
        count+=1
    return docs

def stem_doc_dict(docs):
    c=re.compile("[\[\]\(\)\]!,.=\-;']")
    stemmer =   porter.PorterStemmer()
    engDict =   enchant.Dict("en_US")
    alphan = re.compile("^[a-zA-Z0-9]+")

    for k,x in docs.iteritems():
        x   = re.sub(c,' ',x)
        x2  = ""
        ws  =   re.split("\s+",x)
        for w in ws:
            if not w == '':
                try:
                    if not re.match(alphan,w):
                        pass
                    else:
                        if engDict.check(w):
                            w2  =   stemmer.stem(w)
                            x2  += " "+w2
                        else:
                            x2  += " "+w
                except UnicodeDecodeError:
                    pass
        docs[k] = x2

    return docs
    
    