#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <errno.h>
#include <string.h>
#include <stdbool.h>

static const char * file_name = "intersect_unique.tsv.gz";
unsigned int nrow = 27108124;
int ncol = 121;
int k = 30; 

#define LENGTH 0x2

int dist(int row1, int row2, bool **data) { // Manhattan distance between row1 and row2.
 int i,sum=0;
 for(i=0;i<ncol;i++) {
  sum+= data[row1][i] ^ data[row2][i];
 }
 return(sum);
}

int get_max_i(int *md) { // Returns the index of the element with the maximum value.
 int i, mi, mde = -1;
 for(i=0;i<k;i++) {
  if(md[i] == 0) return(i); // We do NOT expect 0, since we've filtered duplicates. If 0, replace this value.
  if(md[i] > mde) {
   mde = md[i];
   mi = i;
  }
 }
 return(mi);
}

int main() {
 // Allocate memory for the large matrix. It will fit as a bool.
 int i,j;
 bool **data = (bool **)malloc(nrow * sizeof(bool*));//[ncol];// [27108124][121];
 for(i=0;i<nrow;i++) data[i] = (bool*)malloc(ncol*sizeof(bool));;

 // Read in the gzip file.
 gzFile uniq_file;
 uniq_file = gzopen (file_name, "r");

 int bytes_read;
 unsigned char buffer[LENGTH];

 for(i=0;i<nrow;i++) {
  for(j=0;j<ncol;j++) {
   bytes_read = gzread (uniq_file, buffer, LENGTH); // Read two bytes: [0 || 1] and [\t || \n].
   buffer[1] = '\0'; // Set whitespace to string terminator.
   data[i][j] = atoi(buffer); // Not best practice, but seems to work.
   //printf ("%s, %d\n", buffer, data[i][j]);
  }
 }
// printf("Read file.\n");

 // Calculate the minimum distance between rows, and return the row with the minimum distance.
 int d;
 int md[k], mi[k];
 int maxi; // keeps track of the max value we've recorded in md and mi so far.

 for(j=0;j<nrow;j++) {
// j = nrow-1; {
  for(i=0;i<k;i++) {// Go through first k to set md and mi.
   d = dist(j, i, data);
   md[i] = d;
   mi[i] = i;
  }
 
  // Set maxi.
  maxi = get_max_i(md);
  if(md[maxi] == 0) md[maxi] = ncol; // Necessary b/c md[] might be 0.
//  printf("md= %d\tmaxi= %d\n", md[maxi], maxi);
 
  for(i=k;i<nrow;i++) { // Go through remainder of the dataset.
   if(i == j) continue;
   d = dist(j, i, data);
   if(d < md[maxi]) {
    md[maxi] = d;
    mi[maxi] = i;
    
    maxi = get_max_i(md);
   }
  }
//  printf("md= %d\tmaxi= %d\n", md[maxi], maxi);
  for(i=0;i<k;i++) {
   printf("%d\t%d\t%d\n", j, mi[i], md[i]);
  }
//  printf("md= %d\tmaxi= %d\n", md[maxi], maxi);
 }

 gzclose (uniq_file);
 for(i=0;i<nrow;i++) free(data[i]); free(data);
 return(0);

}

