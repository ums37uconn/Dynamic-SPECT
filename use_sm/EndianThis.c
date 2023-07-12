#include <stdio.h>
#include <stdlib.h>


/* ----- Change the Endian format of a file, from one to the other ----- */

int main(int argc, char **argv)
{
  FILE *fidin,*fidout;
  long length=1,num=0,index;
  char data[16];
  
  /* Check arguments */

  if(argc != 4)
  {
    fprintf(stderr, "\nUsage:\n > EndianThis <N> FileToRead.raw FileToWrite.raw\n");
	fprintf(stderr, "\nWhere N is the number of bytes in each data set, should be 2,4 or 8\n");
	exit(1);
  }

  switch(argv[1][0])
  {
    case '2':  length = 2;  break;
    case '4':  length = 4;  break;
    case '8':  length = 8;  break;
   
    default:  
      fprintf(stderr, "The number of bytes in each data set, should be 2,4 or 8, and not \"%s\". \n", argv[2]);
      exit(1);
  }  
 
  if((fidin = fopen(argv[2], "rb"))==NULL)
  { fprintf(stderr, "Could not open file \"%s\" for reading\n",argv[2]);
    exit(1); }

  if((fidout = fopen(argv[3], "wb"))==NULL)
  { fprintf(stderr, "Could not open file \"%s\" for reading\n",argv[3]);
    exit(1); }

  /* read and write files */

  while(fread(&data, sizeof(char), length, fidin)==length)
  {
	num++;
	for(index=length-1;index>=0;index--)
	{
		fputc((int) data[index],fidout);
	}
  }
   
  printf("Wrote %d data sets of length %d.\n",num,length);
  fclose(fidin);
  fclose(fidout);

}
