#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define 	MAX_DIMENSION		100000000
#define 	DEFAULT_BUFFER		100
#define 	MAX_FACTOR		4
#define 	DIMENSIONS_NUM		3
#define 	MAX_INPUT_SYSMATS 	10


void exit_main(int argc, char **argv)
{
  fprintf(stderr, "\nUsage of Matrix Manipulator:\n%s operation (paramaters)\n",argv[0]);
  fprintf(stderr, "The following operations are supported: T for Transpose, S for Shrink, C for Concatenate.\n");
  exit(-1);
}


int main(int argc, char **argv)
{
  if(argc < 2)
  	  exit_main(argc,argv);
  
  switch (argv[1][0]) {
	  case 'T':		return Transpose(argc,argv);
	  case 'S':		return Shrink(argc,argv);
	  case 'C':		return Concat(argc,argv);
	  default:		exit_main(argc,argv);
  	}
}


int Transpose(int argc, char **argv)
{
  int rowsize,colsize;   
  int index,j,mainRowIndex;
  int bufferSize,activeColNum;
  int tmp,tmp2;
  float tmpfl;
  
  FILE *fin,*fout;
  int *Itemp,*outIndex;
  int *currentPos,*activeCellNum,*whichRow;
  int **newI;
  float **newSM;
  long  *rowPointers;
  size_t size_int, size_float, size_long, size_intptr, size_floatptr;
  
  size_int = sizeof(int);
  size_long = sizeof(long);
  size_float = sizeof(float);
  size_intptr = sizeof(int *);
  size_floatptr = sizeof(float *);
	  
  if((argc < 5) || (argc > 6))
  {
    fprintf(stderr, "\nUsage:\n%s T oldSparseMat.file newSparseMat.file <format> (row buffer size)\n",argv[0]);
    fprintf(stderr, "Where format value of 0 is ray driven (column number before row number) and format value of 1 is voxel driven (other way).\n",argv[0]);
    exit(1);
  }
  
  if(argc==6)
  {
	bufferSize=atoi(argv[5]);
	if(bufferSize<=0)
	{  fprintf(stderr, "Error. Buffer size %d must be a positive integer.\n",argv[3]);
	exit(1); }
  }  
  else
	  bufferSize=DEFAULT_BUFFER;
  
  /* Open SparseMat files */
  if((fin = fopen(argv[2], "rb"))==NULL)
   {  fprintf(stderr, "Error. Could not open old sparse matrix file \"%s\" for reading.\n",argv[1]);
      exit(1); }
  
  if(atoi(argv[4])==0)
      {
	      fread(&colsize, size_int, 1, fin);
	      fread(&rowsize, size_int, 1, fin);
      }
  else if(atoi(argv[4])==1)
      {	
	      fread(&rowsize, size_int, 1, fin);
	      fread(&colsize, size_int, 1, fin);	     
      }
  else
      { fprintf(stderr, "Error. Format value %d must be beteen 0 or 1.\n", atoi(argv[4]));
      exit(1);    }
	      
  if((colsize > MAX_DIMENSION) || (colsize <= 0))
  { fprintf(stderr, "Error. Column size in old sparse matrix %d must be beteen 1 and %d.\n", colsize,MAX_DIMENSION);
    exit(1);    }
  if((rowsize > MAX_DIMENSION) || (rowsize <= 0))
  { fprintf(stderr, "Error. Row size in old sparse matrix %d must be beteen 1 and %d.\n", colsize,MAX_DIMENSION);
    exit(1);    }
  
  if((fout = fopen(argv[3], "wb"))==NULL)
  {  fprintf(stderr, "Error. Could not open new sparse matrix file \"%s\" for writing.\n",argv[2]);
     exit(1); }
  fwrite(&rowsize, size_int, 1, fout);
  fwrite(&colsize, size_int, 1, fout);   
 
  rowPointers = (long *)   malloc(size_long  *  colsize);
  
  currentPos = (int *)   malloc(size_int  *  colsize);
  activeCellNum = (int *)   malloc(size_int  *  colsize);
  whichRow = (int *)   malloc(size_int  *  colsize);
  
  Itemp = (int *)   malloc(size_int  *  rowsize);
  outIndex = (int *)   malloc(size_int  *  bufferSize);  

  
  newI = (int **)   malloc(size_intptr  *  bufferSize);
  newSM = (float **) malloc(size_floatptr * bufferSize);

  for(index=0;index<bufferSize;index++)
  {
	 newI[index] = (int *)   malloc(size_int  *  colsize);
	 newSM[index] = (float *) malloc(size_float * colsize); 
  }

  tmpfl=0;
  for(index=0,activeColNum=0;index<colsize;index++)
  {
	  
	  fread(&tmp, size_int, 1, fin);	  
	  if(tmp!=index)
	  { fprintf(stderr, "Error. Old sparse matrix file is corrupt. Expected column %d while reading %d.\n", index,tmp);
	    exit(1);    }
	  fread(&tmp, size_int, 1, fin);
	  
	 
	  if(tmp)
	  {
		  currentPos[activeColNum]=0;
		  activeCellNum[activeColNum]=tmp;
		  whichRow[activeColNum]=index;
		  
		  rowPointers[activeColNum]=ftell(fin);
		  
		  tmpfl+=activeCellNum[activeColNum];
		  fread(&tmp2, size_int, 1, fin);
		  if(tmp2<0)
			  { fprintf(stderr, "Error. First index in column %d has value of: %d.\n", index,tmp2);
			    exit(1);    }
		  for(j=1;j<activeCellNum[activeColNum];j++)
		  {
			fread(&tmp, size_int, 1, fin);
			if(tmp<=tmp2)
			  { fprintf(stderr, "Error. index %d in column %d has value of: %d, after pervious index with value: %d.\n",j, index,tmp,tmp2);
			    exit(1);    }
			tmp2=tmp;
		  }
		  if(tmp>rowsize)
			  { fprintf(stderr, "Error. Last index in column %d with value of: %d is bigger than row size.\n", index,tmp,rowsize);
			    exit(1);    }
		  
		  if(fseek(fin, rowPointers[activeColNum]+8*activeCellNum[activeColNum], SEEK_SET))
		  { fprintf(stderr, "Error. Old sparse matrix file is corrupt. Did not find next column %d.\n", index+1);
		  exit(1);    }
		  activeColNum++;
	  }
  }
  tmpfl /= activeColNum;
  fprintf(stderr, "Read through file successfully. Found %d active columns out of %d, and average use of %.2f rows out of %d. Starting Transpose...\n",activeColNum,colsize,tmpfl,rowsize);
   
  for(mainRowIndex=0;mainRowIndex<rowsize;mainRowIndex+=bufferSize)
  {
	  for(index=0;index<bufferSize;index++)
		  outIndex[index]=0;
	  
	  for(j=0;j<activeColNum;j++)
	  {
		fseek(fin, rowPointers[j]+currentPos[j]*4, SEEK_SET);
		fread(Itemp, size_int, activeCellNum[j]-currentPos[j], fin);
		index=0;		
		fseek(fin, rowPointers[j]+(currentPos[j]+activeCellNum[j])*4, SEEK_SET);
		
		while((Itemp[index]<mainRowIndex+bufferSize) && (currentPos[j]<activeCellNum[j]))
		{
			fread(&tmpfl, size_float, 1, fin);			
			tmp=Itemp[index]-mainRowIndex;
			newI[tmp][outIndex[tmp]]=whichRow[j];
			newSM[tmp][outIndex[tmp]]=tmpfl;
			outIndex[tmp]=outIndex[tmp]+1;
			index++;
			currentPos[j]++;
		}
	  }

	  for(index=0;index<bufferSize;index++)
		if(mainRowIndex+index<rowsize)
	  	{			
			tmp=mainRowIndex+index;
			fwrite(&tmp,size_int, 1, fout);
			fwrite(&outIndex[index], size_int, 1, fout);	
			if(outIndex[index])
			{
				fwrite(newI[index], size_int, outIndex[index], fout);										
				fwrite(newSM[index], size_float, outIndex[index], fout);				
			}			
		}
	fprintf(stderr, ".");
  }
  
  for(index=0;index<bufferSize;index++)
  {
	 free(newI[index]);
	 free(newSM[index]);
  }
  
  free(rowPointers);
  free(currentPos);
  free(activeCellNum);
  free(whichRow);

  free(outIndex);
  free(Itemp);

  free(newI);
  free(newSM);

  fclose(fin);
  fclose(fout);
  
  return 0;
}


int Shrink(int argc, char **argv)
{
  int rowsize,colsize;
  int colDim[DIMENSIONS_NUM],rowDim[DIMENSIONS_NUM];
  int newColDim[DIMENSIONS_NUM],newRowDim[DIMENSIONS_NUM];
  int colFactor[DIMENSIONS_NUM],rowFactor[DIMENSIONS_NUM];
  int colJump[DIMENSIONS_NUM*2],colJump2[DIMENSIONS_NUM];
  int oldRowJump,newRowJump;
  int totColFactor=1,totRowFactor=1;
  int totColDim=1,totRowDim=1;
  int newTotColDim=1,newTotRowDim=1;
  long currentPos=0;
 
  int activeCellNum,activeBasis;
  int dataRead=0,dataWritten=0;
  int index,i1,j1,k1,i2,j2,k2,i3,j3,k3;
  int tmp,tmp2,innerFileIndex,fileIndex;
  float tmpfl;
  FILE *fin,*fout;
  int *Itemp,*hashRow;
  float *SMtemp,*newRow;
  long  *rowPointers;
  size_t size_int, size_float , size_long;
  
  size_int = sizeof(int);
  size_long = sizeof(long);
  size_float = sizeof(float);
  
  
  if(argc != 3+4*DIMENSIONS_NUM+1)
  {
    fprintf(stderr, "\nUsage:\n%s S oldSparseMat.file newSparseMat.file <size of dimension> <factor to divide dimension> .\n",argv[0]);
    fprintf(stderr, "It is assumed there are %d dimensions making up each of the columns and rows (AxBxC=colsize , DxExF=rowsize).\n",DIMENSIONS_NUM);
    exit(1);
  }
  
  for(index=0;index<DIMENSIONS_NUM;index++)
  {
	  colDim[index]=atoi(argv[4+index*2]);
	  rowDim[index]=atoi(argv[4+2*DIMENSIONS_NUM+index*2]);
	  colFactor[index]=atoi(argv[4+index*2+1]);
	  rowFactor[index]=atoi(argv[4+2*DIMENSIONS_NUM+index*2+1]);
	  if((colDim[index] <= 0)||(colDim[index] > MAX_DIMENSION))
	  { fprintf(stderr, "Error. Column dimension (num %d) has value of %d, but must be between 1 and %d.\n",index+1,colDim[index],MAX_DIMENSION);
	    exit(1);    }
	  if((rowDim[index] <= 0)||(rowDim[index] > MAX_DIMENSION))
	  { fprintf(stderr, "Error. Row dimension (num %d) has value of %d, but must be between 1 and %d.\n",index+1,rowDim[index],MAX_DIMENSION);
	    exit(1);    }
	  if((colFactor[index] <= 0)||(colFactor[index] > MAX_FACTOR))
	  { fprintf(stderr, "Error. Column factor for division (num %d) has value of %d, but must be between 1 and %d.\n",index+1,colFactor[index],MAX_FACTOR);
	    exit(1);    }
	  if((rowFactor[index] <= 0)||(rowFactor[index] > MAX_FACTOR))
	  { fprintf(stderr, "Error. Row factor for division (num %d) has value of %d, but must be between 1 and %d.\n",index+1,rowFactor[index],MAX_FACTOR);
	    exit(1);    }
	  newColDim[index]=(int)(colDim[index]/colFactor[index]);
	  newRowDim[index]=(int)(rowDim[index]/rowFactor[index]);
	  totColFactor*=colFactor[index];
	  totRowFactor*=rowFactor[index];
	  totColDim*=colDim[index];
	  totRowDim*=rowDim[index];
	  newTotColDim*=newColDim[index];
	  newTotRowDim*=newRowDim[index];
	  if(newColDim[index]*colFactor[index]!=colDim[index])
		  fprintf(stderr, "Warning. Using column dimension (num %d) with value of %d and factor %d, will create loss of data.\n",index+1,colDim[index],colFactor[index]);
	  if(newRowDim[index]*rowFactor[index]!=rowDim[index])
		  fprintf(stderr, "Warning. Using row dimension (num %d) with value of %d and factor %d, will create loss of data.\n",index+1,rowDim[index],rowFactor[index]);	    
  }
  
  /* Open SparseMat files */
  if((fin = fopen(argv[2], "rb"))==NULL)
   {  fprintf(stderr, "Error. Could not open old SparseMat file \"%s\" for reading.\n",argv[2]);
      exit(1); }
  
  fread(&colsize, size_int, 1, fin);
  fread(&rowsize, size_int, 1, fin);
  
  if(totColDim != colsize)
  { fprintf(stderr, "Error. Multiplied dimensions for column %d do not match colum size from sparse matrix %d.\n", totColDim,colsize);
    exit(1);    }                                                        
  if(totRowDim != rowsize)
  { fprintf(stderr, "Error. Multiplied dimensions for row %d do not match row size from sparse matrix %d.\n", totRowDim,rowsize);
    exit(1);    }
  
   if((fout = fopen(argv[3], "wb"))==NULL)
   {  fprintf(stderr, "Error. Could not open new SparseMat file \"%s\" for writing.\n",argv[3]);
      exit(1); }
   fwrite(&newTotColDim, size_int, 1, fout);
   fwrite(&newTotRowDim, size_int, 1, fout);   
   
   fprintf(stderr, "Using %d columns with %d rows , into %d columns with %d rows .\n", totColDim,totRowDim,newTotColDim,newTotRowDim); 
  
  
  rowPointers = (long *)   malloc(size_long  *  colsize);
  hashRow = (int *) malloc(size_int * rowsize);
  newRow = (float *) malloc(size_float * newTotRowDim);
  Itemp = (int *)   malloc(size_int  *  rowsize);
  SMtemp = (float *) malloc(size_float * rowsize);
  
  
  for(index=0;index<colsize;index++)
  {
	  rowPointers[index]=ftell(fin);
	  fread(&tmp, size_int, 1, fin);
	  if(tmp!=index)
	  { fprintf(stderr, "Error. Sparse matrix file is corrupt. Expected row %d while reading %d. %d.\n", index,tmp,rowPointers[index]);
	    exit(1);    }
	  fread(&tmp, size_int, 1, fin);
	  if(fseek(fin, rowPointers[index]+8+8*tmp, SEEK_SET))
	  { fprintf(stderr, "Error. Sparse matrix file is corrupt. Did not find next row %d.\n", index+1);
	    exit(1);    }
  }
  
  for(i1=0;i1<newRowDim[0];i1++)
	  for(j1=0;j1<newRowDim[1];j1++)
		  for(k1=0;k1<newRowDim[2];k1++)
  		  {
			  newRowJump = i1*newRowDim[1]*newRowDim[2] + j1*newRowDim[2] + k1;
			  for(i2=0;i2<rowFactor[0];i2++)
				  for(j2=0;j2<rowFactor[1];j2++)
					  for(k2=0;k2<rowFactor[2];k2++)
					  {
						  oldRowJump = i1*rowFactor[0]*rowDim[1]*rowDim[2] + j1*rowFactor[1]*rowDim[2] + k1*rowFactor[2]  +  i2*rowDim[1]*rowDim[2] + j2*rowDim[2] + k2;						  	
						  hashRow[oldRowJump] = newRowJump;											  
					  }
		  }
  
  for(i1=0;i1<newColDim[0];i1++)  
  {  
	  colJump[0] = i1*colFactor[0]*colDim[1]*colDim[2]; 
	  colJump2[0]= i1*newColDim[1]*newColDim[2]; 
	  for(j1=0;j1<newColDim[1];j1++)  
	  {
		  colJump[1] = j1*colFactor[1]*colDim[2] + colJump[0];
		  colJump2[1] = j1*newColDim[2] + colJump2[0];
		  for(k1=0;k1<newColDim[2];k1++)  
		  {
			  colJump[2] = k1*colFactor[2] + colJump[1];
			  colJump2[2] = k1+ colJump2[1];
			  for(index=0;index<newTotRowDim;index++) 
				  newRow[index]=0;
			  
			  for(i2=0;i2<colFactor[0];i2++)  
			  {  
				  colJump[3] = i2*colDim[1]*colDim[2] + colJump[2];
				  for(j2=0;j2<colFactor[1];j2++)  
				  {  
					  colJump[4] = j2*colDim[2] + colJump[3];
					  for(k2=0;k2<colFactor[2];k2++)  
					  {  
						colJump[5] = k2 + colJump[4];  					
						fseek(fin, rowPointers[colJump[5]], SEEK_SET);
						
						fread(&tmp, size_int, 1, fin);
						fread(&activeCellNum, size_int, 1, fin);
						
						fread(Itemp, size_int, activeCellNum, fin);
						fread(SMtemp, size_int, activeCellNum, fin);
						
						for(index=0;index<activeCellNum;index++)
							newRow[hashRow[Itemp[index]]]+=SMtemp[index];																	
					  }
				  }
			  }
			  
			  activeCellNum=0;
			  for(index=0;index<newTotRowDim;index++)
				  if(newRow[index]>0)
			  	  {
					Itemp[activeCellNum]=index;
					SMtemp[activeCellNum]=newRow[index]/(totColFactor*totRowFactor);										
					activeCellNum++;
				  }
			  fwrite(&colJump2[2], size_int, 1, fout);
			  fwrite(&activeCellNum, size_int, 1, fout);
			  
			  fwrite(Itemp, size_int, activeCellNum, fout);									
			  fwrite(SMtemp, size_float, activeCellNum, fout);			
		  }
	  }
  }
  free(rowPointers);
  free(hashRow);
  free(newRow);
  free(Itemp);
  free(SMtemp);
  
  fclose(fin);
  fclose(fout);
  
  return 0;
}



int Concat(int argc, char **argv)
{
  int rowsize,colsize;
  int SparseMat_num;
  int activeCellNum;
  int i,j,k;
  int tmp,tmp2,fileIndex;
  float tmpfl;
  FILE *fin,*fout;
  char oldSparseMat_fname[MAX_INPUT_SYSMATS][128], basis_fname[128], newSparseMat_fname[128];
  int *Itemp;
  float *SMtemp;
  size_t size_int, size_float;
  
  size_int = sizeof(int);
  size_float = sizeof(float);
  
  if(argc <= 2)
  {
    fprintf(stderr, "\nUsage:\n%s C #files oldSparseMat.files ... ... newSparseMat.file \n",argv[0]);
    fprintf(stderr, "Where #files is the number of old sparse matrix files.\n",argv[0]);
    exit(1);
  }
  
  SparseMat_num=atoi(argv[2]);
  if((SparseMat_num <= 0)||(SparseMat_num > MAX_INPUT_SYSMATS))
    { fprintf(stderr, "Error. Can not use %d old sparse matrix files, value must be between 1 and %d.\n",SparseMat_num,MAX_INPUT_SYSMATS);
      exit(1);    }
  
  if(argc != 4+SparseMat_num)
  {
    fprintf(stderr, "\nUsage:\n%s #files oldSparseMat.files ... ... newSparseMat.file \n",argv[0]);
    fprintf(stderr, "Where #files is the number of old sparse matrix files.\n",argv[0]);
    exit(1);
  }
  
  for(fileIndex=0;fileIndex<SparseMat_num;fileIndex++)
	  strcpy(&oldSparseMat_fname[fileIndex][0], argv[3+fileIndex]);
  strcpy(newSparseMat_fname, argv[3+SparseMat_num]);
  
  /* Open SparseMat files */
  if((fin = fopen(&oldSparseMat_fname[0][0], "rb"))==NULL)
   {  fprintf(stderr, "Error. Could not open old SparseMat file \"%s\" for reading.\n", &oldSparseMat_fname[0][0]);
      exit(1); }
  if((fout = fopen(newSparseMat_fname, "wb"))==NULL)
   {  fprintf(stderr, "Error. Could not open new SparseMat file \"%s\" for writing.\n", newSparseMat_fname);
      exit(1); }
  
  fread(&colsize, size_int, 1, fin);
  fread(&rowsize, size_int, 1, fin);
  
  fprintf(stderr, "Using %d rows with %d columns.\n", rowsize,colsize);
  
  tmp = colsize*SparseMat_num;
  fwrite(&tmp, size_int, 1, fout);
  tmp = rowsize;
  fwrite(&tmp, size_int, 1, fout);
  
  Itemp = (int *)   malloc(size_int  *  rowsize);
  SMtemp = (float *) malloc(size_float * rowsize);
  
  for(fileIndex=0;fileIndex<SparseMat_num;fileIndex++)
  {
		for(j=0;j<colsize;j++)
		{
			fread(&tmp, size_int, 1, fin);
			fread(&activeCellNum, size_int, 1, fin);
			
			tmp = colsize*fileIndex + j;
			fwrite(&tmp, size_int, 1, fout);
			tmp = activeCellNum;
			fwrite(&tmp, size_int, 1, fout);
						
			fread(Itemp, size_int, activeCellNum, fin);
			fwrite(Itemp, size_int, activeCellNum, fout);
			
			fread(SMtemp, size_float, activeCellNum, fin);
			fwrite(SMtemp, size_float, activeCellNum, fout);
		}
	if(fileIndex<SparseMat_num-1)
	{
		fclose(fin);
		if((fin = fopen(&oldSparseMat_fname[fileIndex+1][0], "rb"))==NULL)
		{  fprintf(stderr, "Error. Could not open old SparseMat file \"%s\" for reading.\n", &oldSparseMat_fname[fileIndex+1][0]);
		exit(1); }
		
		fread(&tmp, size_int, 1, fin);
		fread(&tmp2, size_int, 1, fin);
		if((tmp!=colsize) || (tmp2!=rowsize))
			{  fprintf(stderr, "Error. The dimensions of the newly read SparseMat file %d,%d do not match previous dimensions %d,%d.\n", tmp,tmp2,colsize,rowsize);
			exit(1); }
	}
  }
  free(Itemp);
  free(SMtemp);
  fclose(fin);
  fclose(fout);
  
  return 0;
}

