#include <stdio.h>
#include <stdlib.h>

//#include <metis.h>

//void metisPartMesh_(idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind, idx_t *vwgt, idx_t *vsize, idx_t *nparts, real_t *tpwgts, idx_t *options, idx_t *objval, idx_t *epart, idx_t *npart)
//{
//  int err = 0;
//  
//  err =  METIS_PartMeshNodal(ne, nn, eptr, eind, vwgt, vsize, nparts, tpwgts, options, objval, epart, npart);
//
//  if(err == METIS_OK)
//    printf("Parmetis returns no error!\n");
//
//  if(err == METIS_ERROR)
//    printf("Parmetis returns an error!\n");
//}


void countingsort_(int *array, int*sdata, int *size)
{
  int i, min, max;
  { //calculate min and max
  	min = max = array[0];
	for(i = 1; i < *size; i++) 
	{
	if (array[i] < min)
	min = array[i];
	else if (array[i] > max)
	max = array[i];
	}
  }
  
  int range = max - min + 1;
  int *count = (int *) malloc(range * sizeof(int));	
   if (count == NULL) {
	printf ("Memory allocation error of count.\n");
	return;
  }
   int *sd = (int *) malloc ((*size) * sizeof (int));
   if (sd == NULL) {
	printf ("Memory allocation error of sd.\n");
	return;
  }	
  
  for(i = 0; i < range; i++)
    count[i] = 0;
 
  for(i = 0; i < *size; i++)
    count[ array[i] - min ]++;
  
//   printf ("Range = %d\n", range);
   for (i = 1; i < range; ++i)
 	count [i ] = count [i] + count [i-1];
//   printf ("start executing sd.\n");
   for (i = *size - 1; i > -1; --i) { 
    	sd [count [array[i] - min]-1] = sdata [i];
	count [array[i] - min]--;
   }
//     printf ("data = %d\n", *size);   
   for(i = 0; i < range; i++)
    count[i] = 0;
 
   for(i = 0; i < *size; i++)
   {
       sdata [i] = sd [i];
     count[ array[i] - min ]++;
   }

//   printf ("copying successfully.\n");
  int j, z = 0;
  
  
   for(i = min; i <= max; i++)
    for(j = 0; j < count[ i - min ]; j++)
       array[z++] = i;
//   printf ("copied successfully\n", *size);  
    free(count);
      free (sd);
}

