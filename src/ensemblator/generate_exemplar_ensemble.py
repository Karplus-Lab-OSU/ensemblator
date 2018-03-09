#!/usr/bin/env python2

import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from argparse import ArgumentParser
import numpy as np
from sklearn.cluster import AffinityPropagation
from sklearn import metrics
from scipy.cluster.vq import kmeans, vq
from sklearn.cluster import AgglomerativeClustering
from Bio.PDB import *



# definition for a type of dictionary that will allow new dictionaries to
# be added on demand
class NestedDict(dict):
    def __missing__(self, key):
        self[key] = NestedDict()
        return self[key]


parser = ArgumentParser()
parser.add_argument(
    "-i",
    "--input",
    dest="input",
    type=str,
    metavar="FILE",
    help="This should be a final ensemble prepared by the Ensemblator."
         " Something like 'Global_Overlay_2.5.pdb'"
)
parser.add_argument(
    "-d",
    "--data",
    dest="pairwise",
    type=str,
    metavar="FILE",
    help="This should be a 'pairwise_analysis.tsv' file from an Ensemblator run."
)
parser.add_argument(
    "-o",
    "--output",
    dest="output",
    type=str,
    default="exemplar_ensemble.pdb",
    help="This descriptor will define the final name of the " + "output file."
)
parser.add_argument(
    "--exemplars",
    dest="exemplars",
    type=int,
    default=10,
    help="Select a number of clusters to find, and thus exemplars to include in the final ensemble."
)
parser.add_argument(
    "-c",
    "--color_data",
    dest="color",
    type=str,
    metavar="FILE",
    help="A text file containing one value per row, to use for coloring the t-SNE plots."
)
options = parser.parse_args()



######################


# read in the pairwise data file, we gonna build a matrix
pairwise_file = open(options.pairwise, 'r')

# read every line except the first (the header)
lines = pairwise_file.readlines()[1:]
pairwise_file.close()

# a nested dictionary for each of these three stats we are interested in
dis_score = NestedDict()

# build pairwise dictionaries with all the values of interest
for line in lines:
    x = int(line.split()[0])
    y = int(line.split()[1])

    dis = float(line.split()[5])
    dis_score[x][y] = dis

# this is used to generate ranges, it's so I know what the highest number
# of models is (ie. how many columns and rows I need in my matrix)
max_y = int(max(list(dis_score.keys())))

# making a list that will be formatted in such a way that I can turn it
# into an array or a matrix
dis_score_array_list = []

# go from 0 to the max number of models
for x in range(0, max_y + 1):
    this_x_dis_score = []
    # now do the same for y, now we are going over every pair
    for y in range(0, max_y + 1):

        # fill out the missing values in the array, by duplicating the
        # correct data
        if dis_score[x][y] == {}:
            dis_score[x][y] = dis_score[y][x]

        # append these dictionary values to the list
        this_x_dis_score.append(dis_score[x][y])

    # now append them all, what the list looks like for each x (ie. row)
    # is this [0,23,43,23,53,654,23] where index 0 is x0y0, index 1 is
    # x0y1 etc.
    dis_score_array_list.append(this_x_dis_score)

#### generate co-ocur matrix

co_ocur = NestedDict()

# go from 0 to the max number of models
for x in range(0, max_y + 1):
    # now do the same for y, now we are going over every pair
    for y in range(0, max_y + 1):
        co_ocur[x][y] = 0

print ("\nPerforming iterative Affinity Propagation clustering:\n")

# affinity propagation using dis score
X = np.array(dis_score_array_list)
# convert from a distance matrix to a similarity matrix
X = X * -1

# start with minimum pref
pref = np.min(X[np.nonzero(X)])

af = AffinityPropagation(preference=pref,
                         affinity="precomputed",
                         max_iter=2000).fit(X)
labels = af.labels_
pref = pref * 1.01

# now loop without recording until there are less than n clusters
while len(np.unique(labels)) == len(labels):
    af = AffinityPropagation(preference=pref,
                             affinity="precomputed",
                             max_iter=2000).fit(X)

    labels = af.labels_
    pref = pref * 1.01

# now record and loop until one cluster
while len(np.unique(labels)) != 1:

    # update co-ocur matrix
    # go from 0 to the max number of models
    for x in range(0, max_y + 1):
        # now do the same for y, now we are going over every pair
        for y in range(0, max_y + 1):

            if labels[x] == labels[y]:
                co_ocur[x][y] = co_ocur[x][y] + 1

    af = AffinityPropagation(preference=pref,
                             affinity="precomputed",
                             max_iter=2000).fit(X)

    labels = af.labels_
    pref = pref * 1.01

################## k-means now

print ("\nPerforming iterative K-means clustering:\n")

# get matrixes of these bad boys!
dis_score_asmatrix = np.asmatrix(dis_score_array_list)

# loop from 2 to n-minus-one/2 clusters, some number i times each
for k in range(2, len(labels)):

    i = 0

    while i < 10:

        codebook, dis_score_distortion = kmeans(dis_score_asmatrix, k)
        labels, dist = vq(dis_score_asmatrix, codebook)

        i += 1

        # update co-ocur matrix
        # go from 0 to the max number of models
        for x in range(0, max_y + 1):
            # now do the same for y, now we are going over every pair
            for y in range(0, max_y + 1):

                if labels[x] == labels[y]:
                    co_ocur[x][y] = co_ocur[x][y] + 1

print ("\nPerforming Agglomerative clustering of ensemble of clusters:\n")

# now cluster the co_ocur matrix
co_ocur_array_list = []
# go from 0 to the max number of models
for x in range(0, max_y + 1):
    this_x_co_ocur = []
    # now do the same for y, now we are going over every pair
    for y in range(0, max_y + 1):
        # append these dictionary values to the list
        this_x_co_ocur.append(co_ocur[x][y])

    # now append them all, what the list looks like for each x (ie. row)
    # is this [0,23,43,23,53,654,23] where index 0 is x0y0, index 1 is
    # x0y1 etc.
    co_ocur_array_list.append(this_x_co_ocur)

X = np.array(co_ocur_array_list)

# get the max number of clusters to search for from the user
# optional. Default is 6 (as declared in the option.parser at the top)
k = options.exemplars

labels = AgglomerativeClustering(k, affinity="euclidean", linkage="complete").fit_predict(X)
sil_score = metrics.silhouette_score(X, labels, metric='euclidean')

sil_scores = metrics.silhouette_samples(X, labels, metric='euclidean')
sil_scores_out = open("sil_scores.tsv", "w")
counter = 0
sil_scores_out.write("id" + "\t" +
                     "cluster" + "\t" +
                     "sil_score" + "\n")
groups = {}

for groupID in range(0, len(set(labels))):
    groups[groupID] = []

for label in labels:
    sil_scores_out.write(str(counter) + "\t" +
                         str(label) + "\t" +
                         str(sil_scores[counter]) + "\n")

    groups[label].append(counter)
    counter += 1

# print(groups)


print ("\nThere are " + str(k) + " clusters, with a mean "
       "silhouette score of " + str(sil_score) + ".")
print ("Sillhouette Scores saved in 'sil_scores.tsv'\n")



# now need to calculate the exemplars for each group

mean_dist_dict = NestedDict()

for groupID in groups:
    for modelID in groups[groupID]:
        distances = []
        for y in dis_score[modelID]:
            if y in groups[groupID]:
                distances.append(dis_score[modelID][y])
        mean_dist_dict[groupID][modelID] = np.mean(distances)

exemplar_list = []

for groupID in mean_dist_dict:
    exemplarID = min(mean_dist_dict[groupID], key=mean_dist_dict[groupID].get)
    exemplar_list.append(exemplarID)

print ("Identified models: " + str(exemplar_list) + " as exemplars.")

# now get a pdb file of just the exemplars

pdb_reader = PDBParser(PERMISSIVE=1, QUIET=True)
master_structure = pdb_reader.get_structure("temp", options.input)


class ExemplarSelect(Select):
    def accept_model(self, model):
        if model.get_id() in exemplar_list:
            return 1
        else:
            return 0

io = PDBIO()
io.set_structure(master_structure)
io.save(options.output, ExemplarSelect())

print("Saved final ensemble of exemplar structures as " + str(options.output))

print("Plotting t-SNE plots.")

# init the tsne
tsne = TSNE(n_components=2,
            perplexity=40,
            early_exaggeration=4.0,
            learning_rate=200,
            n_iter=5000,
            n_iter_without_progress=50,
            min_grad_norm=0,
            init='pca',
            method='exact',
            verbose=0)

# run tsne on X, which is the cooccurance matrix generated earlier
reduced = tsne.fit_transform(X)


non_exemplars_x = []
non_exemplars_y = []
exemplars_x = []
exemplars_y = []
exemplars_labels = []
non_exemplars_labels = []

if options.color:
    with open(options.color) as f:
        colors = map(float, f)
    exemplars_colors = []
    non_exemplars_colors = []

title = "exemplar_TSNE"
i = 0
for x in dis_score:
    if x in exemplar_list:
        exemplars_x.append(reduced[i,0])
        exemplars_y.append(reduced[i, 1])
        exemplars_labels.append(labels[i])
        if options.color:
            exemplars_colors.append(colors[i])
    else:
        non_exemplars_x.append(reduced[i, 0])
        non_exemplars_y.append(reduced[i, 1])
        non_exemplars_labels.append(labels[i])
        if options.color:
            non_exemplars_colors.append(colors[i])
    i += 1

plt.figure()
# begin plotting, using the tsne coords in 2d
plt.scatter(non_exemplars_x,
            non_exemplars_y,
            c="blue",
            marker="o",
            edgecolors='none',
            s=80
            )
plt.scatter(exemplars_x,
            exemplars_y,
            c="red",
            marker="D",
            edgecolors='none',
            s=80
            )
# get rid of axis, which is an estimate and should not be taken as gospel
plt.axis('off')
plt.savefig(title + ".svg", bbox_inches='tight')

title = "cluster_TSNE"
plt.figure()
# begin plotting, using the tsne coords in 2d
plt.scatter(non_exemplars_x,
            non_exemplars_y,
            c=non_exemplars_labels,
            marker="o",
            edgecolors='none',
            s=80
            )
plt.scatter(exemplars_x,
            exemplars_y,
            c=exemplars_labels,
            marker="D",
            edgecolors='none',
            s=80
            )
# get rid of axis, which is an estimate and should not be taken as gospel
plt.axis('off')
plt.savefig(title + ".svg", bbox_inches='tight')

if options.color:

    title = "custom_color_TSNE"
    plt.figure()
    plt.scatter(non_exemplars_x,
                non_exemplars_y,
                c=non_exemplars_colors,
                marker="o",
                edgecolors='none',
                s=80,
                cmap=plt.get_cmap('YlOrRd')
                )
    plt.scatter(exemplars_x,
                exemplars_y,
                c=exemplars_colors,
                marker="D",
                edgecolors='none',
                s=80,
                cmap=plt.get_cmap('YlOrRd')
                )

    # get rid of axis, which is an estimate and should not be taken as gospel
    plt.clim(vmin=min(colors), vmax=max(colors))
    plt.colorbar()
    plt.axis('off')
    plt.savefig(title + ".svg", bbox_inches='tight')

print("Done.")


