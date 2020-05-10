import numpy as np
import pickle
from scipy.spatial import cKDTree
from scipy.spatial import distance
from sklearn.neighbors import KDTree, BallTree, NearestNeighbors
import gzip

## Constants
k = 10

# Uniform random distribution
uniform_N = np.random.random((10000, 4))

#mat = np.loadtxt("intersect_unique.tsv.gz", dtype="bool")
mat = np.loadtxt("intersect_all_first_10k.tsv.gz", dtype="bool")

## Jaccard distance.
distance.jaccard(mat[1], mat[1])
distance.jaccard(mat[1], mat[2])

## Fit the kNN. Using Ball Tree algorithm and Jaccard distance.
nn= NearestNeighbors(n_neighbors = k, algorithm='ball_tree', metric='jaccard')
nn.fit(mat)
arr = nn.kneighbors(mat)

## Now write out contents of arr to a file.
out_table = open('knn.intersect_unique.k10.txt', 'w')
#out_table = open('knn.first_10k.k10.txt', 'w')

for x in range(mat.shape[0]):
	for y in range(k):
		silence = out_table.write(str(x) +"\t"+ str(arr[1][x][y]) +"\t"+ str(1 - arr[0][x][y]) +"\n")


out_table.close()

