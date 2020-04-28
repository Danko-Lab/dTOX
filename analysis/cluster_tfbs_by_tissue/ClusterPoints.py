#!/usr/bin/python
import phenograph
import numpy

## Read data.
data = numpy.loadtxt('intersect_all_first_1m.tsv.gz')#, dtype=<class 'integer'>)

## Run phenograph
communities, graph, Q = phenograph.cluster(data)

## Save.
numpy.savetxt('communities.txt.gz', communities, fmt='%.d')
#numpy.savetxt('graph.txt.gz', graph)

