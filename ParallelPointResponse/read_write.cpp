#include<iostream>
#include<cstdlib>
#include<cmath>
#include<fstream>
#include<complex> 
#include <ctime>
#include<string>
#include <sys/stat.h>
#include "cfiMesh2.h"


/*extern "C" {			 			   
		   #include "cfiMesh2.h"
   	    }
*/


using namespace std;

cfiCounter VOXEL_X_COUNT=128;
cfiCounter VOXEL_Y_COUNT=128;
cfiCounter ATN_Z_COUNT=256;

cfiScalar VOXEL_X_WIDTH=4.4;
cfiScalar VOXEL_Y_WIDTH=4.4;
cfiScalar VOXEL_Z_WIDTH=4.4;

cfiCounter MAX_STR_LEN=1024;

cfiDebugStatus g_debugStatus = (cfiProgressDebug|cfiArgDebug|cfiReturnDebug);

cfiCounter *indices, atn_voxel_count;
cfiScalar * origin;

	cfiMesh   *atnMesh;
	cfiArray  *atnMeshDataArray;
	cfiScalar *atnMeshValue;
	cfiMesh   *rayMesh;
	cfiArray  *rayMeshDataArray;
	cfiVector *rayVector;
	cfiScalar *rayValue;
	cfiVector *atnVector;
	cfiScalar *atnVectorValue;


//cfiCounter xAxis=0,yAxis =1, zAxis =2;	


int main(int argc, char **argv){

	FILE *file_ptr;
	char  fileName[MAX_STR_LEN];
	char ATN_FILE_NAME[MAX_STR_LEN];
	char OUT_FILE_NAME[MAX_STR_LEN];

	delete indices;
	indices=new cfiCounter[3];

	indices[xAxis] = VOXEL_X_COUNT;
	indices[yAxis] = VOXEL_Y_COUNT;
	indices[zAxis] = ATN_Z_COUNT;
	atn_voxel_count = VOXEL_X_COUNT * VOXEL_Y_COUNT * ATN_Z_COUNT;

	//cout <<atn_voxel_count<<endl;

	delete origin;
	origin=new cfiScalar[3];
	origin[xAxis] = (-((cfiScalar)(VOXEL_X_COUNT - 1)) / 2.0) * VOXEL_X_WIDTH;
	origin[yAxis] = (-((cfiScalar)(VOXEL_Y_COUNT - 1)) / 2.0) * VOXEL_Y_WIDTH;
	origin[zAxis] = (-((cfiScalar)(ATN_Z_COUNT   - 1)) / 2.0) * VOXEL_Z_WIDTH;

	//cout<<origin[0]<<"\t"<<origin[1]<<"\t"<<origin[2]<<endl;

/*
if ( ( (atnMesh = allocateMesh((cfiMeshType)uniformQuadCornerType, 3, (cfiDataType)scalarType, 3, indices))
		== (cfiMesh *)NULL )
		||   ( (atnMeshDataArray = getMeshDataArrayPtr(atnMesh))
			== (cfiArray *)NULL )
		||   ( (atnMeshValue = (cfiScalar *)getArrayValuePtr(atnMeshDataArray, (cfiCounter *)NULL))
			== (cfiScalar *)NULL ) )
		{
			FAILURE_MSG("cannot allocate attenuation map");
			(void)(exit(cfiFailure));
		}

*/
		strcpy(OUT_FILE_NAME,argv[1]);
		strcpy(ATN_FILE_NAME,argv[3]);

		fprintf(stdout, "reading %s\n", ATN_FILE_NAME);

	if ( ( file_ptr = (FILE *) fopen( ATN_FILE_NAME, "r" ) ) == (FILE *) NULL )
				{
				cout<<"Cannot read attenuation map file: exit(0)"<<endl;
				(void)exit(0);
				}

	//check file size
  struct stat filestatus;
  stat( ATN_FILE_NAME, &filestatus );
  cout<< "total file size:" << filestatus.st_size/1e6 << " Mb\n";
		
	int numberOfATNMapVoxelsRead = (filestatus.st_size)/sizeof( float );
		//= fread( atnMeshValue, sizeof( cfiScalar ), atn_voxel_count, file_ptr );


        if ( numberOfATNMapVoxelsRead != atn_voxel_count ){
		   cout<<"Read number of attenuation map voxels:"<<numberOfATNMapVoxelsRead<<endl;
		   cout<<"Expected number of attenuation map voxels:"<<atn_voxel_count<<endl;	 	
		   cout<<"Mismatch between actual and expected number of attenuation map voxels:exit(1)"<<endl;
		  exit(1);
		}

	cout<<"done!\a"<<endl;
	

}



