#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* ----- Compare two 4-byte float files to get a squared sum of residuals  ----- */

int main(int argc, char **argv)
{
  FILE *file1,*file2,*file3;
  float data1,data2,flag,SSR=0;
  int exitflag = 1,count=0;
  /* Check arguments */

  if(argc != 4)
  {
    fprintf(stderr, "\nUsage:\n > SSR File1.raw File2.raw MaskFile.raw\n");	
    exit(1);
  }

 if((file1 = fopen(argv[1], "rb"))==NULL)
  { fprintf(stderr, "Could not open file \"%s\" for reading\n",argv[1]);
    exit(1); }
  
  if((file2 = fopen(argv[2], "rb"))==NULL)
  { fprintf(stderr, "Could not open file \"%s\" for reading\n",argv[2]);
    exit(1); }

  if((file3 = fopen(argv[3], "rb"))==NULL)
  { fprintf(stderr, "Could not open file \"%s\" for reading\n",argv[3]);
    exit(1); }


  while((fread(&flag, sizeof(float), 1, file3)==1) && (exitflag))
  {
	exitflag = fread(&data1, sizeof(float), 1, file1) * fread(&data2, sizeof(float), 1, file2);	
	if(flag && exitflag)
	{
		SSR += pow(data1-data2,2);
		count++;
	}			
  }
   
  printf("%5.5f\n",pow(SSR/((float)count),0.5));
  fclose(file1);
  fclose(file2);
  fclose(file3);
}

