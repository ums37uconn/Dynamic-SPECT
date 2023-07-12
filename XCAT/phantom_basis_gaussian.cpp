# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
# include <fstream>


using namespace std;

int main(){

int N=360;
int NBasis=360;
float sigma=2;

/*
ofstream test ("respiratory_basis_gaussian_5x5_sigma4.csv");	
		for(int i=0;i<N;i++){
			float x=1.0*i;
			
			for(int j=0;j<NBasis;j++){
			float peak=j;
			float function=exp(-(x-peak)*(x-peak)/sigma);
				test<<function<<",";
				}
			test<<endl;
		}


*/

	double width[N];int count=31;
	srand((unsigned)time(0));
	for(int i=0;i<N;i++){
		//double random_integer = 12.0*rand()*1.0/(RAND_MAX);
		//cout<<i<<"\t" << (int)(random_integer) << endl;
		if(i<=N/2)count--; else count++;
			if(count==0)count=1;
			//cout<<count<<endl;
			width[i]=count;//(int)(random_integer);
		}

ofstream test ("xcat_360basis_sigma2.csv");
 float peak_interval = sigma;
 float peak_width = sigma;
 float peak = sigma;
 
		for(int i=0;i<N;i++){
			float x=1.0*i;
		
			for(int j=0;j<NBasis;j++){
			  peak=j;
		       
			float function=exp(-(x-peak)*(x-peak)/peak_width);
			test<<function<<", ";
				}
			test<<endl;

		}
		cout<<"Number of basis = "<<"\t"<< NBasis<<endl;
		cout<<"gaussian width = "<<"\t"<< sigma <<endl;
		cout<<"Each gaussian has peak interval = "<<"\t"<<peak_interval <<endl;
		

}
