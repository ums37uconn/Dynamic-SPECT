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


//Input/Output file array
float *data;
float *data_new;

string readtype ="4D";//3D or 4D volume data 

FILE *fid;
char sino_fname[128];

int main(int argc, char *argv[]){

 	if (argc != 6 ) {
        std::cout << "Erro:Usage: gaussianFilter4D kernelSize FWHM numberOfTimeFrames inputFilename outputFilename "  << std::endl;
        exit(0);
	}
	int kernelSize = atof(argv[1]); //
    	int FWHM = atof(argv[2]);
    	int numberOfTimeFrames = atof(argv[3]);
	char *inputFilename = argv[4];
	char *outputFilename = argv[5];
	int nTime=numberOfTimeFrames; // number of time frames
	int nVolume4D = nVolume*nTime;
	//int nTime=numberOfTimeFrames; // number of time frames

	//filter specifications
	//int kernelSize =9; //remain odd (5 to 15 reasonable)
	int kernelSizeByTwo = kernelSize/2;
	int kernelSizeCube = kernelSize*kernelSize*kernelSize;
	int kernelSizeSquare = kernelSize*kernelSize; 

	float sigma = FWHM/2.35; 
	float *Kernel;
   


  FILE *fid;
  char data_fname[128];
  strcpy(data_fname,  argv[4]);
  
  data  = new float[nVolume4D];
  data_new  = new float[nVolume4D];
  
  if((fid = fopen(data_fname, "rb"))==NULL)
  {
    fprintf(stderr, "Could not open input data file \"%s\" for reading\n", data_fname);
    exit(1);
  }
  
  //read data file
  fread(data, sizeof(float), nVolume4D, fid);

////////////////////////////
  //Build Kernel
    Kernel = new float[kernelSizeCube];
    int temp = kernelSizeByTwo;
    float s = 2.0 * sigma * sigma; 
    float sum = 0.0;
    for (int i = 0; i <kernelSize;i++) { 
    for (int j = 0; j <kernelSize;j++) { 
    for (int k = 0; k <kernelSize;k++) { 
      float r = sqrt( (i-temp)*(i-temp) + (j-temp)*(j-temp) + (k-temp)*(k-temp) ); 
            Kernel[i*kernelSizeSquare+j*kernelSize+k] = (exp(-(r * r) / s)) / (M_PI * s); 
            sum += Kernel[i*kernelSizeSquare+j*kernelSize+k]; 
      }}} 
    
  //  Kernel Normalization 
  for (int i = 0; i < kernelSize; i++){
    for (int j = 0; j < kernelSize; j++){
      for (int k = 0; k < kernelSize; k++){
          Kernel[i*kernelSizeSquare+j*kernelSize+k] /= sum; 
      }}}
  
 //Activate filter
    int voxelEdgeX = voxelX-kernelSizeByTwo;
    int voxelEdgeY = voxelY-kernelSizeByTwo;
    int voxelEdgeZ = voxelZ-kernelSizeByTwo;
    
  for (int t = 1; t <nTime; t++) { 
                                int s = t*nVolume;
  for (int i = kernelSizeByTwo; i <voxelEdgeX;i++) { 
    for (int j = kernelSizeByTwo; j <voxelEdgeY;j++) { 
      for (int k = kernelSizeByTwo; k <voxelEdgeZ;k++) { 
        
                              int p = (i)*nPlane+(j)*voxelX+(k);
                              int temp = s+p;
        
        for (int l = 0; l<kernelSize; l++) {
          for (int m = 0 ; m<kernelSize; m++) {
            for (int n = 0 ; n<kernelSize; n++) {
              int q = (l)*kernelSizeSquare+(m)*kernelSize+(n);
                data_new[temp] += data[(i+l-kernelSizeByTwo)*nPlane + (j+m-kernelSizeByTwo)*voxelX + (k+n-kernelSizeByTwo) +s] * Kernel[q];
              
          }}}}}}}
  
  FILE * output;
  output = fopen( argv[5], "wb" );
	
  fwrite(data_new, sizeof(float), nVolume4D, output);
 
 tt=clock()-tt;
 //cout<<tt/CLOCKS_PER_SEC<<"seconds"<<endl;
 
  delete[] Kernel;
  delete[]data;
  delete[]data_new;
  
} 

