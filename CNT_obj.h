#ifndef __CNT_OBJ_H_DEFINED__
#define __CNT_OBJ_H_DEFINED__

#include <iostream>
#include <cmath>

#define Pi 3.14159265359

using namespace std;

int seed=-142;

class CNT_obj
{
    private:
		double sig_cnt;
		int N_cnt;
		double *k1_cnt,*k2_cnt,*k3_cnt,*k4_cnt,*tempO_cnt;
		void CNT_ode4(double yold[], double h, double t, void (CNT_obj::*frhs)(double [],double,double[]));
		void CNT_stoerm(double y[], double htot, double ts, int nstep, void (CNT_obj::*frhs)(double [], double, double []));
		void CNT_myrhs(double y[], double t, double f[]);
		void CNT_myrhs_wF(double y[], double t, double f[]);
		double ran1(int &);

	public: 
		int Npoints;
		int i;
		double X0_cnt;
		double d_cnt;
		double K_cnt;
		double kap_cnt;
		double* y_cnt;
		double gam_cnt;
		double temp_cnt;
		double h_cnt;
		double F_cnt;
		double Fper_cnt;
		const static double sigma=7.767E-25;
		const static double kb=1.38E-29;
		
		void timestep( double );
		void timestep_stoermer( double, int);
		void set_F(double );
		void set_gam(double );
		void set_h(double );
		void set_temp(double );
		void set_d(double );
		void update_sig() ;
		void update_Ks();
		double get_sig();
		CNT_obj(double , double , int);

};//end of CNT_obj definition

#endif





