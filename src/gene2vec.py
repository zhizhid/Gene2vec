import gensim, logging
import os
import random
import datetime

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)

print("start!")

sourceDir = "../data" # source directory of the files

#training file format:
#   TOX4 ZNF146
#   TP53BP2 USP12
#   TP53BP2 YRDC

num_db = 0
files = os.listdir(sourceDir)
size = len(files)
gene_pairs = list()
random.shuffle(files)

#load all the data
for fname in files:
    num_db = num_db + 1
    now = datetime.datetime.now()
    print(now)
    print("current file "+ fname + " num: " + str(num_db) + " total files " + str(size))
    f = open(os.path.join(sourceDir, fname), 'r', encoding='windows-1252')
    for line in f:
        gene_pair = line.strip().split()
        gene_pairs.append(gene_pair)
    f.close()

current_time = datetime.datetime.now()
print(current_time)
print("shuffle start " + str(len(gene_pairs)))
random.shuffle(gene_pairs)
current_time = datetime.datetime.now()
print(current_time)
print("shuffle done " + str(len(gene_pairs)))

####training parameters########
dimension = 100 #dimension of the embedding
num_workers = 32 #number of worker threads
sg = 1 #sg =1, skip-gram, sg =0, CBOW
max_iter = 10 #number of iterations
window_size = 1 #The maximum distance between the gene and predicted gene within a gene list

export_dir = "../emb/"

for current_iter in range(1,max_iter+1):
    if current_iter == 1:
        print("gene2vec dimension "+ str(dimension) +" iteration "+ str(current_iter)+ " start")
        model = gensim.models.Word2Vec(gene_pairs, size=dimension, window=window_size, min_count=5, workers=32, iter=1, sg=1)
        model.save(export_dir+"gene2vec_dim_"+str(dimension)+"_iter_"+str(current_iter))
        print("gene2vec dimension "+ str(dimension) +" iteration "+ str(current_iter)+ " done")
        del model
    else:
        current_time = datetime.datetime.now()
        print(current_time)
        print("shuffle start " + str(len(gene_pairs)))
        random.shuffle(gene_pairs)
        current_time = datetime.datetime.now()
        print(current_time)
        print("shuffle done " + str(len(gene_pairs)))

        print("gene2vec dimension " + str(dimension) + " iteration " + str(current_iter) + " start")
        model = gensim.models.Word2Vec.load(export_dir+"gene2vec_dim_"+str(dimension)+"_iter_"+str(current_iter-1))
        model.train(gene_pairs,total_examples=model.corpus_count,epochs=model.iter)
        model.save(export_dir+"gene2vec_dim_"+str(dimension)+"_iter_"+str(current_iter))
        print("gene2vec dimension " + str(dimension) + " iteration " + str(current_iter) + " start")
        del model