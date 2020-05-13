#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <errno.h>
#include <string.h>
#include <stdbool.h>
#include "knnd.h"
#include "vec.h"

static const char * file_name = "intersect_unique_sort.tsv.gz";
static const char * out_file_name = "knnd_table.tsv.gz";
/*unsigned*/ int nrow = 1000;// 27108124; // Test on a smaller dataset.
int ncol = 121;
int k = 10;

#define LENGTH 0x2

float jaccard(bool* row1, bool* row2, int d) { // Manhattan distance between row1 and row2.

 int and=0,or=0; // MUCH Faster to keep these conditions separate.
 for(int i=0;i<d;i++) {
  and+= row1[i] && row2[i];
  or += row1[i] || row2[i];
 }

 return( 1.0f - (float)and/(float)or ); // If 3/4 match in two nodes w/ 4
}

int read_data(dataset_t* data) {
 data->values = malloc((sizeof(bool*) * nrow) + (nrow * ncol * sizeof(bool)));
 data->size = nrow;
 data->dim = ncol;

 printf("Reading the gzip file.\n");
 gzFile uniq_file;
 uniq_file = gzopen (file_name, "r");

 int bytes_read;
 unsigned char buffer[LENGTH];

 for(int i=0;i<nrow;i++) {
  data->values[i] = (bool*) (data->values + data->size) + i * data->dim;
  for(int j=0;j<ncol;j++) {
   bytes_read = gzread (uniq_file, buffer, LENGTH); // Read two bytes: [0 || 1] and [\t || \n].
   buffer[1] = '\0'; // Set whitespace to string terminator.
   data->values[i][j] = atof(buffer); // Not best practice, but seems to work.
  }
 }

 gzclose(uniq_file);
 return(0);
}

int main() {
 // Allocate memory for the large matrix. It will fit as a bool.
 printf("setting up variables.\n");

 // Read in the gzip file.
 dataset_t data;
 if (read_data(&data)) return 1;

 // Calculate the minimum distance between rows, and return the row with the minimum distance.
 vec_t* B = nn_descent(data, &jaccard, k, 1.0, 0.001);

 // Print out the heap.
 gzFile out_file;
 out_file = gzopen (out_file_name, "w");

 for(int i=0;i<data.size;i++) {
  gzprintf(out_file, "%d\t%d\t%f\n", i, i, 0.0);

  for(int j=0;j<B[i].size;j++) {
   gzprintf(out_file, "%d\t%d\t%f\n", i, B[i].arr[j].id, B[i].arr[j].val);
  }
 }
 gzclose(out_file);
 
 heap_list_free(B, data.size);
 free(data.values);
 return(0);
}

