#!/usr/bin/python
import sys
import gzip

in_file = gzip.open('intersect_final.bed.gz')

out_unique = gzip.open('intersect_unique.tsv.gz', 'wb')
out_IDgroup= gzip.open('intersect_identity_group.tsv.gz', 'wb') # This is 0-base line # in intersect_unique.

## New dictionary.
count = 0
uniq_lines = dict()

## Loop through the file.
for line in in_file:
	# Check if line has already been observed.
	if line in uniq_lines:
		out_IDgroup.write(uniq_lines[line]+"\n")
	else: 
		out_IDgroup.write(str(count)+"\n")
		out_unique.write(line)
		uniq_lines[line] = str(count)
		count+=1

in_file.close()

out_unique.close()
out_IDgroup.close()
