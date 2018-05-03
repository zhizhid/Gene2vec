import sys
import codecs
import numpy as np
import matplotlib.pyplot as plt
import plotly.plotly as py
import plotly.graph_objs as go

from sklearn.manifold import TSNE

Y = np.loadtxt("TSNE_tsne_data_m_i_100_total.txt")

rdVB = np.loadtxt("TSNE_label_m_i_100_total.txt", dtype='str')

geneEmbedSet = set()

with open("GSE_mat/GSE_200_mat", 'r') as file:
    for line in file:
        geneEmbedSet.add((line.split("\t")[0]))
file.close()
data = list()
with open("MsigDB/simMSigTop10.txt", 'r') as readFile:
    for line in readFile:
        rdVBList = list(rdVB.tolist())
        geneList = list()
        tmpList = line.split("\t")
        name = tmpList[0]
        indexList = list()
        n=len(tmpList)
        for i in range(2,n-1):
            if tmpList[i] in geneEmbedSet:
                geneList.append(tmpList[i])
        for str in geneList:
            if str in rdVBList:
                indexList.append(rdVBList.index(str))
        trace = go.Scatter(
            x=Y[indexList, 0],
            y=Y[indexList, 1],
            #    z = Y[:, 2],
            mode='markers',
            text=rdVB[indexList],
            hoverinfo="text",
            name = name)
        cpTrace = trace
        data.append(cpTrace)

        # update tsne index
        np.delete(Y, indexList)
        np.delete(rdVB, indexList)
        # refresh varaiables
        del trace
        tmpList.clear()
        geneList.clear()
        indexList.clear()
readFile.close()

trace = go.Scatter(
            x=Y[:2000, 0],
            y=Y[:2000, 1],
            #    z = Y[:, 2],
            mode='markers',
            text=rdVB[:2000],
            hoverinfo="text",
            marker = dict(
                color = 'lightgrey',)
        )
cpTrace = trace
tmpTrace = data[0]
data[0] = cpTrace
data.append(tmpTrace)
# plotly.offline.iplot({
#     "data": [Scatter(x=reduced_matrix[:, 0], y=reduced_matrix[:, 1])],
#     "layout": Layout(title="hello world")
# })

layout = go.Layout(
    showlegend=True
)
fig = go.Figure(data=data, layout=layout)
plot_url = py.plot(fig, filename='tsne')