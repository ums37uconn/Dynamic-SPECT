#include<iostream>
#include<cstdlib>
#include<cmath>
#include<fstream>
#include<ctime>
#include<cstring>

using namespace std;

float tt = clock();
//input image specifications

int voxelX=128;//number of voxel in X
int voxelY=128;// in Y
int voxelZ=128;// in Z
int nVolume = voxelX*voxelY*voxelZ;
int nPlane = voxelX*voxelY;
int sumStartFrame =60;// start summing from sumStartFrame

float totalActivity = 10*3.7e1; //in MBq
float scanTime = 1200;// total scan time (sec)
float voxelSize = 0.44*0.44*0.44;



float *data;
float *data_new;
float *data_sum;

string readtype ="4D";

FILE *fid;
char sino_fname[128];

int main(int argc, char *argv[]){

 	if (argc != 5 ) {
        std::cout << "Error:Usage: sum4D numberOfTimeFrames inputFilename outputFilename "  << std::endl;
        exit(0);
	}
    	int numberOfTimeFrames = atof(argv[1]);
	char *inputFilename = argv[2];
	char *outputFilename = argv[3];
	char *outputSumFilename = argv[4];
	int nTime=numberOfTimeFrames; //number of time frames
	int nVolume4D = nVolume*nTime;
	int timeFramesDiff = numberOfTimeFrames - sumStartFrame;

//read input data file
  FILE *fid;
  char data_fname[128];
  strcpy(data_fname,  argv[2]);
  
  data  = new float[nVolume4D];
  data_new  = new float[nVolume4D];
  data_sum  = new float[nVolume];
  
  if((fid = fopen(data_fname, "rb"))==NULL)
  {
    fprintf(stderr, "Could not open input data file \"%s\" for reading\n", data_fname);
    exit(1);
  }

  fread(data, sizeof(float), nVolume4D, fid);

float totalCount = 0;
      
  	for (int i = 0; i <nVolume;i++) {
		for (int t = sumStartFrame; t <nTime; t++) {      
			totalCount+= data[t*nVolume + i];  
          }}
	
	float countRate = totalCount/scanTime;
	float calibrationFactor = (countRate/totalActivity)*voxelSize;
/*
	cout<<totalActivity<<"\t"
		<<scanTime<<"\t"
		<<totalCount<<"\t"
		<<countRate<<"\t"
		<<calibrationFactor<<"\t"
		<<timeFramesDiff<<"\t"
		<<endl;
*/
//normalization of 4D data
for (int i = 0; i <nVolume;i++) {
		for (int t = 0; t <nTime; t++) {      
			data_new[t*nVolume + i]=data[t*nVolume + i]/calibrationFactor;   
          }}

  FILE * output;
  output = fopen( argv[3], "wb" );	
  fwrite(data_new, sizeof(float), nVolume4D, output);



//sum from frame sumStartFrame to endFrame
for (int i = 0; i <nVolume;i++) {
		for (int t = sumStartFrame; t <nTime; t++) {     
                	data_sum[i] += data_new[t*nVolume + i]; 
          }}

//time average
for (int i = 0; i <nVolume;i++) {   
                	data_sum[i] /= timeFramesDiff; 
          }

  FILE * outputSum;
  outputSum = fopen( argv[4], "wb" );
	
  fwrite(data_sum, sizeof(float), nVolume, outputSum);
 

/*
int count=0;
for (int i = 0; i <nVolume;i++) {     
			if (data_new[i] > 1){ 
				//cout<<  data_new[i] <<endl; 
				count++;
          }}
  
cout<<count*1.0/nVolume<<endl;
*/
 
 //tt=clock()-tt;
//cout<<tt/CLOCKS_PER_SEC<<"seconds"<<endl;
 

  delete[]data;
  delete[]data_new;
  delete[]data_sum;
  
} 

