#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define MAX_BASIS_NUM 5000
#define MAX_BASIS_LEN 50000
#define MAX_INPUT_SYSMATS 10

int main(int argc, char **argv)
{

  float *basis[MAX_BASIS_NUM];
  int basisNum,basisLen;
  int rowsize,colsize;
  int single_proj_size;
  int sysmat_num;
  int activeCellNum,activeBasis;
  int dataRead=0,dataWritten=0;
  int i,j,k,l;
  int tmp,tmp2,innerFileIndex,fileIndex;
  float tmpfl;
  FILE *fin,*fout;
  char oldsysmat_fname[MAX_INPUT_SYSMATS][128], basis_fname[128], newsysmat_fname[128];
  int *Itemp;
  float *SMtemp;
  size_t size_int, size_float;
  
  size_int = sizeof(int);
  size_float = sizeof(float);
  
  if(argc <= 1)
  {
    fprintf(stderr, "\nUsage:\n > %s <#1> oldsysmat.files ... ... basis.file newsysmat.file <#2> \n",argv[0]);
    fprintf(stderr, " > Where #1 is the number of old system matrix files, and #2 is the size of a single projection (one angle). \n",argv[0]);
    exit(1);
  }
  
  sysmat_num=atoi(argv[1]);
    if((sysmat_num <= 0)||(sysmat_num > MAX_INPUT_SYSMATS))
    { fprintf(stderr, "Error. Can not use %d old system matrix files, value must be between 1 and %d.\n",sysmat_num,MAX_INPUT_SYSMATS);
      exit(1);    }
  
  if(argc != 5+sysmat_num)
  {
    fprintf(stderr, "\nUsage:\n > %s <#1> oldsysmat.files ... ... basis.file newsysmat.file <#2> \n",argv[0]);
    fprintf(stderr, " > Where #1 is the number of old system matrix files, and #2 is the size of a single projection (one angle). \n",argv[0]);
    exit(1);
  }
  
  
  for(fileIndex=0;fileIndex<sysmat_num;fileIndex++)
	  strcpy(&oldsysmat_fname[fileIndex][0], argv[2+fileIndex]);
  strcpy(basis_fname,   argv[2+sysmat_num]);
  strcpy(newsysmat_fname, argv[3+sysmat_num]);
  
  
  // -- read basis file  
  if((fin = fopen(basis_fname, "rb"))==NULL)
   { fprintf(stderr, "Error. Could not open basis file \"%s\" for reading.\n", basis_fname);
    exit(1);    }
  dataRead+=fread(&basisNum, size_int, 1, fin);
  dataRead+=fread(&basisLen, size_int, 1, fin);
  if((basisLen <= 0)||(basisNum <= 0)||(basisLen > MAX_BASIS_LEN)||(basisNum > MAX_BASIS_NUM))
   { fprintf(stderr, "Error. Parameters of Bspline file are out of bounds,  basis num: %d, basis length: %d.\n", basisNum,basisLen);
    exit(1);    }
       
  // -- read basis coefficients
  for(i=0;i<basisNum;i++)
  {
	basis[i] = (float *)   malloc(size_float  *  basisLen);  
  	for(j=0;j<basisLen;j++)	 
		dataRead+=fread(&basis[i][j], size_float, 1, fin);
  }
  fclose(fin);
  if ( dataRead != (basisNum*basisLen+2))	 
   {  fprintf(stderr, "Error. Bspline file corrupt %d != %d.\n",dataRead,(basisNum*basisLen+2));     
     exit(1); 	}

  // open sysmat files
  if((fin = fopen(&oldsysmat_fname[0][0], "rb"))==NULL)
   {  fprintf(stderr, "Error. Could not open old sysmat file \"%s\" for reading.\n", &oldsysmat_fname[0][0]);
      exit(1); }
  if((fout = fopen(newsysmat_fname, "wb"))==NULL)
   {  fprintf(stderr, "Error. Could not open new sysmat file \"%s\" for writing.\n", newsysmat_fname);
      exit(1); }
  
  fread(&colsize, size_int, 1, fin);
  fread(&rowsize, size_int, 1, fin);
  innerFileIndex=0;
  
  single_proj_size=atoi(argv[4+sysmat_num]);
    if((single_proj_size <=0)||(single_proj_size > colsize))
    { fprintf(stderr, "Error. Single projection size %d is not valid for column size %d.\n", single_proj_size,colsize);
      exit(1);    }
  
  fprintf(stderr, "Using %d rows with %d columns, single projection size is: %d.\n", rowsize,colsize,single_proj_size);
  
  tmp=basisLen*single_proj_size*sysmat_num;
  fwrite(&tmp, size_int, 1, fout);
  tmp=basisNum*rowsize;
  fwrite(&tmp, size_int, 1, fout);
  
  Itemp = (int *)   malloc(size_int  *  rowsize);
  SMtemp = (float *) malloc(size_float * rowsize);
 
  for(fileIndex=0;fileIndex<sysmat_num;fileIndex++)
  {
	for(i=0;i<basisLen;i++)
  	{
		dataRead=0;
		dataWritten=0;
		activeBasis=0;
		for(k=0;k<basisNum;k++)
			if(basis[k][i]>0)
				activeBasis++;
		
		for(j=0;j<single_proj_size;j++)
		{
			fread(&tmp, size_int, 1, fin);
			fread(&activeCellNum, size_int, 1, fin);
			
			tmp=basisLen*single_proj_size*fileIndex+single_proj_size*i+j;
			fwrite(&tmp, size_int, 1, fout);
			tmp=activeBasis*activeCellNum;
			fwrite(&tmp, size_int, 1, fout);
			
			if(activeCellNum>rowsize)
			{ fprintf(stderr, "Error. Single projection size %d is not valid for row size %d.\n", single_proj_size,rowsize);
			exit(1);    }
			
			else if(activeCellNum>0)
			{		
				dataRead+=fread(Itemp, size_int, activeCellNum, fin);
				for(k=0;k<basisNum;k++)
					if(basis[k][i]>0)
						for(l=0;l<activeCellNum;l++)
						{
							tmp=Itemp[l]+k*rowsize;
							dataWritten+=fwrite(&tmp, size_int, 1, fout);	
						}
				
				dataRead+=fread(SMtemp, size_int, activeCellNum, fin);
				for(k=0;k<basisNum;k++)
					if(basis[k][i]>0)
						for(l=0;l<activeCellNum;l++)
						{
							tmpfl=SMtemp[l]*basis[k][i];
							dataWritten+=fwrite(&tmpfl, size_float, 1, fout);	
						}
			} 
		}
		
		if(dataRead*activeBasis!=dataWritten)
		{ fprintf(stderr, "Error. Number of read (%d*%d=%d) and written (%d) cells does not match.\n", dataRead,activeBasis,dataRead*activeBasis,dataWritten);
		exit(1);    }
		
		if((innerFileIndex+2*single_proj_size)<=colsize)
			innerFileIndex+=single_proj_size;
		else
		{					
			fclose(fin);
			if((fin = fopen(&oldsysmat_fname[fileIndex][0], "rb"))==NULL)
			{  fprintf(stderr, "Error. Could not open old sysmat file \"%s\" for reading.\n", &oldsysmat_fname[fileIndex][0]);
			exit(1); }
			fread(&tmp, size_int, 1, fin);
			fread(&tmp, size_int, 1, fin);  
			innerFileIndex=0;		
		}
		
	}
	if(fileIndex<sysmat_num-1)
	{
		fclose(fin);
		if((fin = fopen(&oldsysmat_fname[fileIndex+1][0], "rb"))==NULL)
		{  fprintf(stderr, "Error. Could not open old sysmat file \"%s\" for reading.\n", &oldsysmat_fname[fileIndex+1][0]);
		exit(1); }
		 fread(&tmp, size_int, 1, fin);
		 fread(&tmp2, size_int, 1, fin);
		 if((tmp!=colsize)||(tmp2!=rowsize))
			{  fprintf(stderr, "Error. The dimensions of the newly read sysmat file %d,%d do not match previous dimensions %d,%d.\n", tmp,tmp2,colsize,rowsize);
			exit(1); }
		 innerFileIndex=0;
	}
  }
  
  free(Itemp);
  free(SMtemp);
  for(i=0;i<basisNum;i++)  
	free(basis[i]);  
  
  fclose(fin);
  fclose(fout);
  
  return 1;
}
