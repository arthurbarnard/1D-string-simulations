#include "CNT_obj.cpp"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <sstream>
#include <ctime>

using namespace std;

#define Pi 3.14159265359
	
int main(int argc, char**argv)
{	
	/*load in parameters at runtime, a typical call is like the following:
	CNT_run 1 initial_condition_filename output_filename seedTable.txt
	
	This would mean the input filename would be "initial_condition_filename_001.txt"
	and the output filename would be "output_filename_001.txt".
	The random seed would come from the second number in the seedTable, since enumeration
	starts at 0.
	*/
	int fileNum = atoi(argv[1]);
	char *fileIn = argv[2];
	char *fileOut = argv[3];
	char *seedIn = argv[4];
	char *massIn= argv[5];
	char *KIn= argv[6];
	char *kapIn= argv[7];
	
	//instantiate and initialize some parameters
	long long int  n,i,istep, nstep=4000000, nwait,splot=1500, Nin;
	double t,y,tinit=0, T=80000*64, h=T/nstep,d,L,temperature, dispRate, Nmass, N_K_cnt, N_kap_cnt;
	double*yin, *ymass, *y_K_cnt, *y_kap_cnt;
	//myRead will be used to load in binary data files
	double* myRead(char[],long long int&);
	
	/*-------edit simulation parameters here------------*/
	nstep=4.5E9; //total number of time steps
	nwait=10E6; //number of time steps for initial thermalization
	h=1.28; //time step in ps
	d=3; //diameter in nm
	L=10000.; //length in nm
	temperature=9.0; //temperature in Kelvin
	dispRate=100E-7; //dissipation rate in (ps^-1)
	
	T=nstep*h; 
	
    //initialize initfile name and load in data
	char *initfile=new char[20000];
	for(i=0;i<20000;i++)initfile[i]='\0';
	strcat(initfile,fileIn);
	strcat(initfile,"_000");
	int S2=strlen(initfile); 
	initfile[S2-1]=fmod(fileNum,10.0)+48;	
	initfile[S2-2]=fmod(fileNum/10.0,10.0)+48;
	initfile[S2-3]=fmod(fileNum/100.0,10.0)+48;
	strcat(initfile,".txt");
	cout<<initfile<<endl;
	yin=(double*)myRead(initfile,Nin);
	 
    //instantiate CNT_obj with the diameter (d) and free length (L) as inputs
    //(Nin/6) is the number masses in the input file. L needs to be correct since
    //the raw positions in the input file don't tell you the free length of the CNT
	CNT_obj myCNT=CNT_obj(d,L,Nin/6);
	
	//initialize and load in mass, extensional rigidity, and bending rigidity files
	char *massfile=new char[20000];
	for(i=0;i<20000;i++)massfile[i]='\0';
	strcat(massfile,massIn);
	strcat(massfile,"_000");
	int S3=strlen(massfile); 
	massfile[S3-1]=fmod(fileNum,10.0)+48;	
	massfile[S3-2]=fmod(fileNum/10.0,10.0)+48;
	massfile[S3-3]=fmod(fileNum/100.0,10.0)+48;
	strcat(massfile,".txt");
	ymass=(double*)myCNT.massRead(massfile,Nmass);

	char *Kfile=new char[20000];
	for(i=0;i<20000;i++)Kfile[i]='\0';
	strcat(Kfile,KIn);
	strcat(Kfile,"_000");
	int S4=strlen(Kfile); 
	Kfile[S4-1]=fmod(fileNum,10.0)+48;	
	Kfile[S4-2]=fmod(fileNum/10.0,10.0)+48;
	Kfile[S4-3]=fmod(fileNum/100.0,10.0)+48;
	strcat(Kfile,".txt");
	y_K_cnt=(double*)myCNT.K_cnt_read(Kfile,N_K_cnt);
	
	char *kapfile=new char[20000];
	for(i=0;i<20000;i++)kapfile[i]='\0';
	strcat(kapfile,kapIn);
	strcat(kapfile,"_000");
	int S5=strlen(kapfile); 
	kapfile[S5-1]=fmod(fileNum,10.0)+48;	
	kapfile[S5-2]=fmod(fileNum/10.0,10.0)+48;
	kapfile[S5-3]=fmod(fileNum/100.0,10.0)+48;
	strcat(kapfile,".txt");
	y_kap_cnt=(double*)myCNT.kap_cnt_read(kapfile,N_kap_cnt);
	
    int middleInd = floor(Nin/6/2-1)*6; //index pertaining to the middle of the CNT
    myCNT.set_h(h);
	myCNT.set_temp(temperature);
	myCNT.set_gam(100E-7); 
	
	//load the input file data into the CNT_obj
	for(i=0;i<myCNT.Npoints*6;i++)
	{
		myCNT.y_cnt[i]=yin[i];
	}
	
	for(i=0;i<myCNT.Npoints;i++)
	{
		myCNT.massarray[i]=ymass[i];
		myCNT.K_cnt_array[i]=y_K_cnt[i];
		myCNT.kap_cnt_array[i]=y_kap_cnt[i];;
	}
	
	//initialize outputFile name from input arguments
	char *outputFile=new char[20000];
	for(i=0;i<20000;i++)outputFile[i]='\0';
	strcat(outputFile,fileOut);
	strcat(outputFile,"_000");
    int Sl=strlen(outputFile);
   	outputFile[Sl-1]=fmod(fileNum,10.0)+48;	
	outputFile[Sl-2]=fmod(fileNum/10.0,10.0)+48;
	outputFile[Sl-3]=fmod(fileNum/100.0,10.0)+48;
	//open main output file to store full data stream
	strcat(outputFile,".txt");
	ofstream file ((&outputFile[0]), ios::out| ios::binary);
	outputFile[Sl]='\0';
	//open output file to store data on a single point in the middle of the stream
	strcat((&outputFile[0]),"Y1.txt");
	ofstream fileY ((&outputFile[0]), ios::out| ios::binary);
	outputFile[Sl]='\0';
	
	//read in a file with a series of random integer seed for the random number generator
	//the specific seed is picked out of a list based on the fileNum parameter
	ifstream seedFile (seedIn, ios_base::in | ios_base::binary);
	seedFile.seekg(fileNum*4,ios::beg);
	seedFile.read((char*)(&seed),sizeof(int));
	seed*=-1;
		
	//main simulation loop
	tinit=0;
	for(i=0;i<nstep;i++)
	{	
		//set the time of each step
		t=tinit+h*i;
		//turn off dissipation and make a small step in time
		if(i==nwait) 	
		{
			myCNT.set_gam(0); 
			myCNT.set_h(h/5); 
			myCNT.timestep(t);
			myCNT.set_h(h); 
		}
		//periodically store data and print repsentative points to terminal
		if(i%splot==0)
		{
			cout<<myCNT.y_cnt[middleInd]<<"   "
				<<myCNT.y_cnt[middleInd+1]<<"   "
				<<myCNT.y_cnt[middleInd+2]<<endl;
				
			y=0;
			
			for(int k=0;k<myCNT.Npoints*6;k++) 
			{
				y=myCNT.y_cnt[k];	
				//writes data to a binary file in the order of x,y,z,Vx,Vy,Vz for each
				//localized mass
				file.write((char*)(&y),8);
				//writes a single point along a single axis (Y) to a separate file
				//for quick analysis
				if(k==(middleInd+7)) fileY.write((char*)(&y),8);
			}
		}
		//uses runge kutta to prethermalize the system before nwait, 
		//uses stoermer algorithm with no dissipation to simulate hamiltonian dynamics
		if(i>nwait)
			myCNT.timestep_stoermer(t,5);
		else
			myCNT.timestep(t);
 	}
			
	file.close();
	fileY.close();
	
	return 0;

	return EXIT_SUCCESS;
}//end main

// reads in a binary file and creates a double array
double* myRead(char filename[],long long int &n){
	ifstream file (filename, ios::in);
	file.seekg(0,ios::beg);
	double size=(double)(-file.tellg());
	file.seekg(0,ios::end);
	size+=(double)file.tellg();
	char* h=(char*)malloc((int)size);
	file.seekg(0,ios::beg);
	file.read(h,size);
	double* x=(double*)h;
	n=(int)size/8;
	file.close();
	return x;
}












