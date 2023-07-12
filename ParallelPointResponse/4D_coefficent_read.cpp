#include<iostream>
#include<cstdlib>
#include<cmath>
#include<fstream>
#include <ctime>
#include<string>
#include <sys/stat.h>

typedef int cfiCounter;
typedef float cfiScalar;

using namespace std;

int StringLength=128;
cfiScalar *CoefficientValues;

int main(){

	FILE *file_ptr;
	char Coefficient[StringLength];
	strcpy(Coefficient,"/home/uttam/use_sm/recon/coefficient.raw");
	if ( ( file_ptr = (FILE *) fopen( Coefficient, "r" ) ) == (FILE *) NULL ){
				cout<<"Cannot read coefficient file:exited"<<endl;
				exit(1);
				}
	//check file size
 	 struct stat filestatus;
 		stat( Coefficient, &filestatus );
 		cout<< "total file size:" << filestatus.st_size/1e6 << " Mb\n";
	//read file

	CoefficientValues=new cfiScalar[111990];
	cfiScalar count=fread(CoefficientValues, sizeof( cfiScalar ),111990, file_ptr );
	for(int i=0;i<119000;i++)
	cout<<CoefficientValues[i]<<endl;
	delete[]CoefficientValues;
}



