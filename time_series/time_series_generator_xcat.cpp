#include <iostream>
#include<cstdlib>
#include <fstream>
#include <vector>
#include <cmath>


using namespace std;


/*
Henry Wang, original author
Fares Alhassen, editor
Uttam Shrestha, modified extensively
Physics Research Laboratory
UCSF, 2012
*/

/*!
 @param int size
 @param const float
 @return none
 */

/*
void writeFloatToBinary (int size, const std::vector< float > &writeData)
{
    FILE * output;
    output = fopen("rawdata.f32le", "wb");
    fwrite( &writeData, size, sizeof(float), output);
}
*/


//linear interpolation function
float interpolate(float inputPoint, int functionNumber, int number0fSplinePerFunction , const float *splinePoints )
{
    float youtput;
    int xinput = floor(inputPoint);
    int basisAdjusted = xinput + 2;//2 points a begining are not spline points
    //std::cout << splinePoints[4] << std::endl;

    //adjusting?
    int functionPlace = (functionNumber - 1) * number0fSplinePerFunction;
    basisAdjusted = basisAdjusted + functionPlace;

    //based on formula http://i.ajdesigner.com/cdn/linear_interpolation_equation.png
    youtput = splinePoints[basisAdjusted] + ((inputPoint- xinput)*(splinePoints[basisAdjusted+1] - splinePoints[basisAdjusted]));

    return youtput;
}


float getNonuniformTimeRange(float start, float step)
{
    float timeNow = start;
    for(int i = 0; i > 6; i++)
    {
        start = start + step;
    }

}



//start time, endtime, number of times, data file, spline file
int main(int argc, char *argv[])
{
    if (argc != 7 ) {
        std::cout << "Usage: time_series_generator start_projection end_projection number_of_samples data_file spline_file output_file"  << std::endl;
        exit(0);
    }

    float secondsTotal = 0;


    float startTime = atof(argv[1]);
    float endTime = atof(argv[2]);
    int numberOfTimes = atof(argv[3]);

    bool debug = false; //to select slice
    int numberofXVoxels = 128;
    int numberofYVoxels = 128;
    int numberofZVoxels = 128;

   int nVolume =numberofXVoxels*numberofYVoxels*numberofZVoxels;

    int startSlice = 0;
    int endSlice = 0;

    int numberOfSlicesToWrite = endSlice - startSlice + 1;
    int numberOfVoxelsPerSlice = numberofXVoxels * numberofYVoxels;
    int numberOfVoxelsPerSet = numberOfVoxelsPerSlice * numberofZVoxels;

    int numberofVoxelsPerSelectedSet = numberOfVoxelsPerSlice * numberOfSlicesToWrite;

    //std::cout << numberofVoxelsPerSelectedSet  << std::endl;

    int startVoxel = numberOfVoxelsPerSlice * startSlice;
    int endVoxel = startVoxel + numberOfVoxelsPerSlice * numberOfSlicesToWrite;

    std::ifstream::pos_type inputDataSizeInBytes;
    char *memblock;
    char *outputFilename = argv[6];

    //read in binary to float
    //std::ifstream file ("recon_Pt1_str_dynamic_sm_no_pr_128x128x128x6_mh.f32le.iter054", std::ios::in|std::ios::binary|std::ios::ate);
    //std::ifstream file ("impulse_128x128x128x5_4.raw", std::ios::in|std::ios::binary|std::ios::ate);
    //s2 = new string (argv[4]);
    //std::cout << argv[4]  << std::endl;
    std::ifstream file ( argv[4], std::ios::in|std::ios::binary|std::ios::ate);
    if ( file.is_open() ) {
        inputDataSizeInBytes = file.tellg();
        memblock = new char [ inputDataSizeInBytes ];
        file.seekg(0, std::ios::beg);
        file.read( memblock, inputDataSizeInBytes );
        file.close();
        file.clear();
    }
    else {
         printf("Unable to open file\n");
         exit(0);
    }

    float* inputData = reinterpret_cast< float* >( memblock );
    int numberOfdata = (inputDataSizeInBytes / sizeof( float ));

    std::ifstream::pos_type inputSplineSizeInBytes;
    char *memblockSplines;

    //std::ifstream Splinefile ("singlehead_spline_fitted.bin", std::ios::in|std::ios::binary|std::ios::ate);
    //std::ifstream Splinefile ("time_basis_five_splines_432.bin", std::ios::in|std::ios::binary|std::ios::ate);

    std::ifstream Splinefile (argv[5], std::ios::in|std::ios::binary|std::ios::ate);
    if ( Splinefile.is_open() ) {
        inputSplineSizeInBytes = Splinefile.tellg();
        memblockSplines = new char [ inputSplineSizeInBytes ];
        Splinefile.seekg(0, std::ios::beg);
        Splinefile.read( memblockSplines, inputSplineSizeInBytes );
        Splinefile.close();
        Splinefile.clear();
   	 } 
	else
    	printf("Unable to open file\n");

    //casting
    int* basisAndPoint = reinterpret_cast< int* >( memblockSplines );
    float* basisData = reinterpret_cast< float* >( memblockSplines );
    int numberOfSplinedata = (inputSplineSizeInBytes / sizeof( float )) - 2;

    int numberOfBasisFunctions = 0;
    int numberofTimePoints = 0;
    numberOfBasisFunctions = basisAndPoint[0];
    numberofTimePoints = basisAndPoint[1];

cout<<"numberOfBasisFunctions="<<numberOfBasisFunctions<<"\t"
	<<"numberofTimePoints="<<numberofTimePoints<<endl;

    //std::cout << numberOfSplinedata/numberOfBasisFunctions  << std::endl;
    //std::cout << numberofTimePoints  << std::endl;
    //std::cout << numberOfSplinedata  << std::endl;

    int entriesPerCoefficient = numberOfdata/numberOfBasisFunctions;

//cout<<entriesPerCoefficient<<endl;

    //typedef std::vector< std::vector<float> > vectorOfVectorVoxels;
    int sizeX, sizeY;
    int count = 0;

    sizeY = entriesPerCoefficient;
    sizeX = numberOfTimes;

    float	basisStore[numberOfBasisFunctions];
    float   coeffTimesBasis[numberOfBasisFunctions];

    float timeNow = startTime;
    float timePointNow = 0;

    secondsTotal = endTime - startTime;
    float timePointPerSecond = numberofTimePoints/secondsTotal;
    float timeIncrement = (secondsTotal)/numberOfTimes;


//for(int i=0;i<numberofTimePoints*numberOfBasisFunctions;i++)
//cout<<basisData[i]<<endl;
//for(int i=0;i<numberOfdata;i++)
//if(inputData[i]>0.2)cout<<inputData[i]<<endl;


    float sum = 0;

	//writeFloatToBinary(sizeX*sizeY, voxelSeries);
	FILE * output;
	output = fopen( outputFilename, "wb" );

/*
    for ( int x = 0; x < sizeX; x++) {
        timeNow = startTime + (x*timeIncrement);
        timePointNow = timeNow*(timePointPerSecond);
        //std::cout <<  timePointNow  << std::endl;
        for  (int z = 1; z <= numberOfBasisFunctions; z++) {
            		basisStore[z-1]= interpolate(timePointNow, z, numberofTimePoints, basisData);
        	}        
        
        for(int y = 0; y < sizeY; y++) {
            //y+2*entriesPerCoefficient for example adjusts the input data to the correct place
            count = 0;
            sum = 0;
            	for (int a = 1; a <= numberOfBasisFunctions; a++) {
                 	sum += inputData[y+count*entriesPerCoefficient] * basisStore[a-1];
                	//std::cout <<  count  << std::endl;
               	 	count++;
           		 }

*/
//////older version for the whole volume
	 
	for ( int t = startTime; t<endTime; t++) {//time
		if(t%1==0)
 			for ( int i = 0; i<nVolume; i++) {//time
	
			sum = 0;
            		for (int j = 0; j < numberOfBasisFunctions; j++) {
                 		sum += inputData[j*nVolume+i] * basisData[j*numberofTimePoints+t];
           			 }

            		float *sump;

            		sump = &sum;

            	fwrite(sump, sizeof(float), 1, output);


	}}
          
	/*

//new version with X,Y,Z ROI
for ( int t = startTime; t<endTime; t++) {//time
		if(t%5==0)
 			//for ( int i = 0; i<nVolume; i++) {

			for ( int x = numberofXVoxels/4+10; x<3*numberofXVoxels/4; x++) {
				for ( int y = numberofYVoxels/4+10; y<3*numberofYVoxels/4; y++) {
					for ( int z = numberofZVoxels/4+10; z<3*numberofZVoxels/4; z++) {
	
			int i=x*numberofYVoxels*numberofZVoxels+y*numberofZVoxels+z;	

			sum = 0;
            		for (int j = 0; j < numberOfBasisFunctions; j++) {
                 		sum += inputData[j*nVolume+i] * basisData[j*numberofTimePoints+t];
           			 }

            		float *sump;

            		sump = &sum;

            	fwrite(sump, sizeof(float), 1, output);

	}}}//for volume

	}
	*/

    delete[] memblock;
    delete[] memblockSplines;
    //delete[] voxelSeries1;

}

