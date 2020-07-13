#include "CNT_obj.h"
#include <fstream>
#include <cstdlib>

	// CNT_ode4 is an implemenetation of the 4th order runge kutta method. The right hand side (frhs)
	// stores the details of the forces in the equation of motion, this is a separate equation called below.	
	// ---adapted from Numerical Recipes in C---
	void CNT_obj::CNT_ode4(double yold[], double h, double t, void (CNT_obj::*frhs)(double [],double,double[]))
	{ 
	  	int i;
		double t2=t+h*0.5, t3=t+h;
	
		(this->*frhs)(yold,t,k1_cnt);
		for(i=0;i<N_cnt;i++) tempO_cnt[i]=yold[i]+0.5*h*k1_cnt[i];
		(this->*frhs)(tempO_cnt,t2,k2_cnt);
		for(i=0;i<N_cnt;i++) tempO_cnt[i]=yold[i]+0.5*h*k2_cnt[i];
		(this->*frhs)(tempO_cnt,t2,k3_cnt);
		for(i=0;i<N_cnt;i++) tempO_cnt[i]=yold[i]+h*k3_cnt[i];
		(this->*frhs)(tempO_cnt,t3,k4_cnt);
		
		for(i=0;i<N_cnt;i++) 
		{
		 yold[i]+=h*(k1_cnt[i]+2.0*k2_cnt[i]+2.0*k3_cnt[i]+k4_cnt[i])/6.0;
			//	cout<<"inside CNT_ode4   "<<i<<endl;
		}
	     return;
	}//end ode4

	// CNT_stoerm is an implemenetation of Stoermers aglorithm for 2nd order conservative eautions. The right hand side (frhs)
	// stores the details of the forces in the equation of motion, this is a separate equation called below.
	// ---adapted from Numerical Recipes in C---
	void CNT_obj::CNT_stoerm(double y[], double htot, double ts, int nstep, void (CNT_obj::*frhs)(double [], double, double []))
	{
		int i,j,nn,n,neqns;
		double h,h2,halfh,t;
		
		h=htot/nstep;
		(this->*frhs)(y,ts,k1_cnt);
		
		halfh=0.5*h;
		neqns=N_cnt/2;
		for (i=12;i<N_cnt-12;i+=6) {
			for(j=0;j<3;j++){
				n=i+j;
				tempO_cnt[n]=y[n]+(tempO_cnt[n+3]=h*(y[n+3]+halfh*k1_cnt[n+3]));
			}
		}
		t=ts+h;
		(this->*frhs)(tempO_cnt,t,k3_cnt);
		h2=h*h;
		for (nn=1;nn<nstep;nn++) {
			for (i=12;i<N_cnt-12;i+=6){
				for (j=0;j<3;j++){
					n=i+j;
					tempO_cnt[n] += (tempO_cnt[n+3] += h2*k3_cnt[n+3]);
				}
			}
			t += h;
		(this->*frhs)(tempO_cnt,t,k3_cnt);
		}
		for (i=12;i<N_cnt-12;i+=6) {
			for (j=0;j<3;j++){
				n=i+j;
				k3_cnt[n+3]=tempO_cnt[n+3]/h+halfh*k3_cnt[n+3];
				k3_cnt[n]=tempO_cnt[n];
			}
		}
		for(i=12;i<N_cnt-12;i++)
		{
			y[i]=k3_cnt[i];
		}
	}


	/*-------------------------------------------------------*/
	// myrhs is the specific set of 1st order ODEs to be solved
	void CNT_obj::CNT_myrhs(double y[], double t, double f[])
	{
	    double d0,d1,d2,d3,d01,d11,d22,d32,d12,temp,rr,rr1,rr2,d[N_cnt/6],Del[N_cnt],a;
	    int m,i,n;
	    
// 	    double* massRead(char [], long long); When you compile this program, we get an undefined reference error. I believe these defintions are left unused so I commented them out
// 	    double* K_cnt_read(char[], long long);
// 	    double* kap_cnt_read(char[], long long);
	    char* initfile=new char[20000];
		
		
	    //measure some lengths
	    for(n=0;n<N_cnt/6-1;n++)
	    {
	        d[n]=0;
	        for(i=0;i<3;i++)
		{
	            temp=y[(n+1)*6+i]-y[n*6+i];
	            Del[n*6+i]=temp;
	            temp*=temp;
	            d[n]+=temp;
		}
	        d[n]=sqrt(d[n]);
	    }
	    //perform some dot-products
		rr=Del[6]*Del[0]+Del[7]*Del[1]+Del[8]*Del[2];
		rr1=Del[12]*Del[6]+Del[13]*Del[7]+Del[14]*Del[8];
		rr2=Del[12]*Del[18]+Del[13]*Del[19]+Del[14]*Del[20];

		//calculates the first order differentials for each coordinate, to be integrated elsewhere
		for(n=12;n<N_cnt-12;n+=6){
			// dx/dt = Vx ... dz/dt = Vz
			f[n]=y[n+3];
			f[n+1]=y[n+4];
			f[n+2]=y[n+5];	
			
			//pick relevant distances and and pre-multiply pairs
	        d0=d[n/6-2];
	        d1=d[n/6-1];
	        d01=d1*d0;
	        d11=d1*d1;
	        d2=d[n/6];
	        d12=d1*d2;
	        d22=d2*d2;
	        d3=d[n/6+1];
	        d32=d3*d2;
	        // dVx/dt = -gamma*Vx + F_Vx(t) ... dVz/dt = -gamma*Vz + F_Vz(t)
	        
	    for(i=0;i<3;i++){
	            m=n+i;
	            f[m+3]=0*(-gam_cnt*y[m+3])/(massarray[n/6])+ //damping forces
			     (K_cnt_array[n/6]/massarray[n/6])*(Del[m]*(1-X0_cnt/d2)-Del[m-6]*(1-X0_cnt/d1))+ //forces due to stretching
		    	     (kap_cnt_array[n/6-1]/massarray[n/6])*(Del[m-12]-Del[m-6]*rr/d11)/d01-
			     (kap_cnt_array[n/6+1]/massarray[n/6])*(Del[m+6]-Del[m]*rr2/d22)/d32+
			     (kap_cnt_array[n/6]/massarray[n/6])*(Del[m]*(1+rr1/d22)-Del[m-6]*(1+rr1/d11))/d12;//forces due to bending
			
			}
	        //shift stored values for reuse on next data point
		rr=rr1;
	        rr1=rr2;
	        rr2=Del[n+6]*Del[n+12]+Del[n+7]*Del[n+13]+Del[n+8]*Del[n+14];			
	    }
	}
	//same as above, but with a constant force applied in the neg. z-direction
	void CNT_obj::CNT_myrhs_wF(double y[], double t, double f[])
	{
		CNT_myrhs(y, t, f);
		for(int i=12;i<N_cnt-12;i+=6) f[i+5]-=Fper_cnt;
	}
	//random number generator adapted from numerical recipes in C
	double CNT_obj::ran1(int &idum)
	{
		const int IA=16807, IM=2147483647, IQ=127773, IR=2836, NTAB=32;
		const int NDIV=(1+(IM-1)/NTAB);
		const double EPS=3.0e-16,AM=1.0/IM,RNMX=(1.0-EPS);
		static int iy=0;
		static int iv[NTAB];
		int j,k;
		double temp;
		
		if (idum <=0 || !iy){
			if (-idum <1) idum=1;
			else idum = -idum;
			for(j=NTAB+7;j>=0;j--){
				k=idum/IQ;
				idum=IA*(idum-k*IQ)-IR*k;
				if (idum<0) idum += IM;
				if (j<NTAB) iv[j] = idum;
			}
			iy=iv[0];
		}
		k=idum/IQ;
		idum=IA*(idum-k*IQ)-IR*k;
		if (idum<0) idum +=IM;
		j=iy/NDIV;
		iy=iv[j];
		iv[j]=idum;
		if ((temp=AM*iy)>RNMX) return RNMX;
		else return temp;
	}

	//performs one step of numerical integration with runge kutta
	void CNT_obj::timestep( double t_in)
	{	
		//implements thermal noise as velocity perturbations
		for(i=15;i<N_cnt-12;i+=6)
            {
                //y_cnt[i]+=sig_cnt*(2*ran1(seed)-1);
                //y_cnt[i+1]+=sig_cnt*(2*ran1(seed)-1);
                //y_cnt[i+2]+=sig_cnt*(2*ran1(seed)-1);
			}
		//numerical inegration performed	
		(this->CNT_ode4)(y_cnt, h_cnt, t_in, &CNT_obj::CNT_myrhs_wF);
	
	}

	//performs one step of numerical integration with Stoermers method
	void CNT_obj::timestep_stoermer( double t_in,int nstep)
	{
			 (this->CNT_stoerm)(y_cnt, h_cnt, t_in, nstep , &CNT_obj::CNT_myrhs_wF);	
	}

	//update a constant force applied
	void CNT_obj::set_F(double F_in)
	{
		F_cnt=F_in;
		Fper_cnt=F_in/((double)Npoints-3.0);
	}

	//changes the dissipation rate
	void CNT_obj::set_gam(double gam_in)
	{
		gam_cnt=gam_in;
		update_sig();
	}

	//changes the discretized unit of time
	void CNT_obj::set_h(double h_in)
	{
		h_cnt=h_in;
		update_sig();
	}

	//sets the temperature
	void CNT_obj::set_temp(double temp_in)
	{
		temp_cnt=temp_in;
		update_sig();
	}

	//changes the CNT diameter
	void CNT_obj::set_d(double d_in)
	{
		 d_cnt=d_in;
		 update_Ks();
		 update_sig();
	}
	
	//changes "velocity perturbation amplitude" based on the simulation parameters
	void CNT_obj::update_sig() 
	{
		sig_cnt=sqrt(h_cnt*temp_cnt/sigma*kb/Pi/d_cnt/X0_cnt)*sqrt(gam_cnt*6.0);
	}

	//update spring constants based on hollow CNT model 
	void CNT_obj::update_Ks()
	{
		
		K_cnt=446/X0_cnt/X0_cnt;
		kap_cnt=K_cnt*d_cnt*d_cnt/8.0;
	}

	//returns the "velocity perturbation amplitude"
	double CNT_obj::get_sig()
 	{
  	 	   return sig_cnt;
 	}
 	
 	void CNT_obj::load_Property_arrays()
 	{
 		
 		cout<<N_cnt<<endl;
 		ymass=(double*)massRead("", Nmass);
		for(int n=0;n<N_cnt/6;n++)
		{
			massarray[n]=ymass[n];
		}
		
		y_K_cnt=(double*)K_cnt_read("", N_K_cnt);
		for (int n=0;n<N_cnt/6;n++)
		{
			K_cnt_array[n]=y_K_cnt[n];
		
		}
		
		y_kap_cnt=(double*)kap_cnt_read("", N_kap_cnt);
		for (int n=0;n<N_cnt/6;n++)
		{
			kap_cnt_array[n]=y_kap_cnt[n];
			cout<<n<<"  "<<massarray[n]<<"  "<<K_cnt_array[n]<<"  "<<kap_cnt_array[n]<<endl;
		}
	 }
 	
	/*-------------------------------------------------------*/
	//class definition of CNT_obj that stores relevant simulation parameters and coordinates of internal particsls
	
	CNT_obj::CNT_obj(double d_in, double L_in, int Npoints_in)
	{
	    Npoints=Npoints_in;
		N_cnt=Npoints_in*6;
		k1_cnt=new double[5*N_cnt];
		k2_cnt=k1_cnt+N_cnt;
		k3_cnt=k2_cnt+N_cnt;
		k4_cnt=k3_cnt+N_cnt;
		tempO_cnt=k4_cnt+N_cnt;
		
		cout<<"test in CNT_obj"<<endl;
		X0_cnt=L_in/((double)Npoints_in-3.0);   

		set_d(d_in);
		y_cnt=new double[N_cnt];
		massarray=new double[Npoints_in];
		K_cnt_array=new double[Npoints_in];
		kap_cnt_array=new double[Npoints_in];
		
		for(int i=0;i<N_cnt;i++) y_cnt[i]=0;
		for(int i=0;i<N_cnt;i+=6) y_cnt[i]=i/6*X0_cnt;
		void CNT_myrhs(double [], double, double []);
		void CNT_ode4(double [], double, double, void(double [], double, double[]));
		void CNT_stoerm(double [], double, double, int, void(double [], double, double[]));
		
		load_Property_arrays();
		cout<<N_cnt<<endl;
		cout<<massarray[0];
		

	}

double* CNT_obj::massRead(char filename[],long long Nmass)
{
	ifstream massfile("mass.txt", ios::in);
	massfile.seekg(0,ios::beg);
	double sizemass=double(-massfile.tellg());
	massfile.seekg(0,ios::end);
	sizemass+=double(massfile.tellg());
	char* h=(char*)malloc((int)sizemass);
	massfile.seekg(0,ios::beg);
	massfile.read(h,sizemass);
	double* x=(double*)h;
	Nmass=(int)sizemass/8;
	massfile.close();
	return x;	
}

double* CNT_obj::K_cnt_read(char filename[], long long N_K_cnt)
{
	ifstream K_cnt_file("K_cnt.txt",ios::in);
	K_cnt_file.seekg(0,ios::beg);
	double size_K_cnt=double(-K_cnt_file.tellg());
	K_cnt_file.seekg(0,ios::end);
	size_K_cnt+=double(K_cnt_file.tellg());
	char* h=(char*)malloc((int)size_K_cnt);
	K_cnt_file.seekg(0,ios::beg);
	K_cnt_file.read(h,size_K_cnt);
	double* x=(double*)h;
	N_K_cnt=(int)size_K_cnt/8;
	K_cnt_file.close();
	return x;
}

double* CNT_obj::kap_cnt_read(char filename[], long long N_kap_cnt)
{
	ifstream kap_cnt_file("kap_cnt.txt",ios::in);
	kap_cnt_file.seekg(0,ios::beg);
	double size_kap_cnt=double(-kap_cnt_file.tellg());
	kap_cnt_file.seekg(0,ios::end);
	size_kap_cnt+=double(kap_cnt_file.tellg());
	char* h=(char*)malloc((int)size_kap_cnt);
	kap_cnt_file.seekg(0,ios::beg);
	kap_cnt_file.read(h,size_kap_cnt);
	double* x=(double*)h;
	N_kap_cnt=(int)size_kap_cnt/8;
	kap_cnt_file.close();
	return x;
}

	










