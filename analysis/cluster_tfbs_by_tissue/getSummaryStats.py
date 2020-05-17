#!/usr/bin/python
import sys
import gzip
from array import *

in_file = gzip.open('intersect_final.bed.gz')

## New dictionary.
rowCounts = [0] * 121
colCounts = [0] * 121

## Loop through the file.
for line in in_file:
	value = csv.reader(fd, delimiter="\t")
	

	# Check if line has already been observed.
	if line in uniq_lines:
		out_IDgroup.write(uniq_lines[line]+"\n")
	else: 
		out_IDgroup.write(str(count)+"\n")
		out_unique.write(line)

in_file.close()

rowCountsFile = gzip.open('rowCounts.intersect_final.txt.gz', 'wb')
colCountsFile = gzip.open('colCounts.intersect_final.txt.gz', 'wb')



rowCountsFile.close()
colCountsFile.close()


