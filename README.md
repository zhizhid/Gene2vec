# Gene2vec
Codes repository for paper "Gene2Vec: Distributed Representation of Genes Based on Co-Expression".
This work will be presented at [ICIBM](http://icibm2018.zhaobioinfo.org/).

We've provided a pre-trained gene2vec file that achieved best clusteredness in our experiment, see pre_trained_emb/gene2vec_dim_200_iter_9.txt.
## Installing Gene2vec

### Requirements
Gene2vec relies on Python 3.6, TensorFlow 1.6+, gensim 3.4.0, numpy 1.14.0+, and matplotlib 2.1.2+

The Multicore t-SNE relies on a [GitHub repository](https://github.com/DmitryUlyanov/Multicore-TSNE).
Please follow their installation instruction.

### Install
To download the codes, please do:

```
git clone https://github.com/jingcheng-du/Gene2vec.git
cd Gene2vec/
pip install -r requirements.txt
```

## Usage

### gene2vec.py
This script takes gene co-expression data as the input and output binary gene2vec file. It relies on gensim Word2vec module.
The detailed parameters for Word2vec module can be seen in their [document](https://radimrehurek.com/gensim/models/word2vec.html).

For co-expression data files, each line should have a pair of genes that is separated by a space.
The format of the files should be like:
```
TOX4 ZNF146
TP53BP2 USP12
TP53BP2 YRDC
```
As the output is in binary format, here we provide the generateMatrix.py that can
convert the binary embedding to txt format. Each line starts with the gene name and is followed by the vector.

You can run this code directly. However, to speed up the training, please make sure you have a C compiler before installing gensim, see this
[article](https://rare-technologies.com/word2vec-in-python-part-two-optimizing/).

### tsne_multi.py
We provide this script to map the high dimensional gene2vec to 2-D array.
It relies on on a Parallel t-SNE implementation [GitHub repository](https://github.com/DmitryUlyanov/Multicore-TSNE).

It takes the gene2vec txt file as the input and generates two separate files:
a label file that each line is the gene name and a data file that each line has 2-D vector for the corresponding gene.

We implement a mutlipe-thread tsne training here. You can specify a set of iterations you'd like to train in parallel:
```
p = Pool(6)
p.map(TSNEWoker, ["100","5000","10000","20000","50000","100000"])
```

### GGIPNN_Classification.py
GGIPNN stands for gene-gene interaction predictor neural network. Given training, validation and testing dataset for gene-gene interaction,
this script will leverage multilayer perceptron (defined in GGIPNN.py) to predict gene-gene interaction and generate the AUC score.

Each dataset (training, validation and testing) should have two files:

1. a text file that each line has the names of the gene pair:
```
GPNMB BAP1
GPR34 CARD16
ELF5 TGFB2
LILRB2 NCR2
```

2. a label file that each line has either 1 (the corresponding gene pair has the gene-gene interaction) or 0 (doesn't have the gene-gene interaction):
```
0
0
1
1
```

You can specify whether to initialize the embedding layer using pre-trained gene2vec or initialize the embedding layer randomly.
You can also specify whether the embedding layer is trainable during the training process. Just to choose True or False here:
```
tf.flags.DEFINE_boolean("use_pre_trained_gene2vec", True, "use_pre_trained_gene2vec")
tf.flags.DEFINE_boolean("train_embedding", False, "train_embedding")
```

### GTRxFigure.py
This script takes gene2vec tsne files (label and data) and tissue-specific genes expression files as the input
and generates tissues-specific genes experssion maps.

The tissue-specific genes expression file should follow the following format, each line starts with the gene name and followed by the z score:
```
ATRX	0.598962527411
TCOF1	-0.264690317639
NSRP1	0.716336803551
OPA1	0.223145913678
RHEB	0.978021549909
SEMA5B	-0.110593590242
```

## Citation

If you use these codes in your publications, please cite this [paper](https://www.biorxiv.org/content/early/2018/03/22/286096):

```
@article{2018gene2vec,
  title={{Gene2Vec}: Distributed Representation of Genes Based on Co-Expression},
  author={Du, Jingcheng and Jia, Peilin and Dai, YuLin and Tao, Cui and Zhao, Zhongming and Zhi, Degui},
  journal={bioRxiv},
  year={2018}
}
```
