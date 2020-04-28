#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <errno.h>
#include <string.h>
#include <stdbool.h>

static const char * file_name = "intersect_unique.tsv.gz";
unsigned int nrow = 27108124;
int ncol = 121;
int k = 60;

#define LENGTH 0x2

inline int dist(int row1, int row2, bool **data) { // Manhattan distance between row1 and row2.
 int i,xor=0;
 for(i=0;i<ncol;i++) {
  xor+= data[row1][i] ^  data[row2][i];
 }
 if(xor >= 5) return -1;

 int and=0,or=0; // MUCH Faster to keep these conditions separate.
 for(i=0;i<ncol;i++) {
  and+= data[row1][i] && data[row2][i];
  or += data[row1][i] || data[row2][i];
 }

// if(row1 == row2 && (float)and/(float)or < 0.6)
//  printf("row1: %d\trow2: %d\txor: %d\tand: %d\tor: %d\t and/or: %f\n", row1, row2, xor, and, or, (float)and/(float)or);
 return( ((float)and/(float)or >= 0.6)?xor:-1 ); // If 3/4 match in two nodes w/ 4
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
  }
 }

 // Calculate the minimum distance between rows, and return the row with the minimum distance.
 int d, count;

 for(j=0;j<nrow;j++) {
  count = 1;
  for(i=j;i<nrow;i++) {// Go through j .. nrow
   d = dist(j, i, data);
   if(d>=0) { // d>0 prevents reporting connections to self. >=0 reports self connections.
     printf("%d\t%d\t%d\n", j, i, d);
     count++;
     if(count > k) break;
   }
  }
 }
 
 gzclose (uniq_file);
 for(i=0;i<nrow;i++) free(data[i]); free(data);
 return(0);

}

