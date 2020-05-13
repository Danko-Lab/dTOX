#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <errno.h>
#include <string.h>
//#include <stdbool.h>
#include "knnd.h"
#include "vec.h"

static const char * file_name = "intersect_unique_sort.tsv.gz";
unsigned int nrow = 10000;// 27108124; // Test on a smaller dataset.
int ncol = 121;
int k = 30;
int maxRows = 1000000;

#define LENGTH 0x2

inline float jaccard(float* row1, float* row2, int d) { // Manhattan distance between row1 and row2.
// float acc = 0.0f, and = 0.0f, or = 0.0f;

 int and=0,or=0; // MUCH Faster to keep these conditions separate.
 for(i=0;i<d;i++) {
  and+= row1[i] && row2[i];
  or += row1[i] || row2[i];
 }

 return( (float)and/(float)or ); // If 3/4 match in two nodes w/ 4
}

int main(int argc, char** argv) {
 // Allocate memory for the large matrix. It will fit as a bool.
 int i,j;
 dataset_t data;
 data.values = malloc((sizeof(float*) * nrow) + (nrow * ncol * sizeof(float))); 
 
 data.size = nrow;
 data.dim = ncol;


 // Read in the gzip file.
 gzFile uniq_file;
 uniq_file = gzopen (file_name, "r");

 int bytes_read;
 unsigned char buffer[LENGTH];

 for(i=0;i<nrow;i++) {
  data->values[i] = (float*) (data->values + ncol) + i * d;
  for(j=0;j<ncol;j++) {
   bytes_read = gzread (uniq_file, buffer, LENGTH); // Read two bytes: [0 || 1] and [\t || \n].
   buffer[1] = '\0'; // Set whitespace to string terminator.
   data.values[i][j] = atof(buffer); // Not best practice, but seems to work.
  }
 }

 // Calculate the minimum distance between rows, and return the row with the minimum distance.
 vec_t* B = nn_descent(data, &jaccard, k, 1.0, 0.001);

 // Print out the heap.
 for(i=0;i<data.size;i++) {
  heap_check(&B[i], 0);
  for(j=0;j<&B[i].size;j++) {
   printf("%d\t%d\t%f\n", i, &B[i].arr[j].id, &B[i].arr[j].val);
  }
 }
 
 gzclose(uniq_file);
 heap_list_free(B, data.size);
 free(data.values);
 return(0);
}

