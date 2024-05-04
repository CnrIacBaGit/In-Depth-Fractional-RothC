// Solver for fractional-order zRothC model
// author - Vsevolod Bohaienko
#include "gamma.h"

//////////////////////////////////////
///// global variables ///////////////
//////////////////////////////////////
const double HtoPmult=9.8;

int N=50;  // number of nodes
int testing_mode=0; // testing with constructed solution (1), C equation without advection (2)
int debug_level=0;
int zoutstep=1;
// linear solver (TFQMR)
int ls_max_iter=30; // maximal number of iterations
double ls_eps=1e-12; // accuracy threshold
double fr_eps=5e-5;
double ls_min_tau=0.001; // minimal time step length, s
double max_tau=1*24*3600.0; // maximal time step length, s
double ls_mult=1.25; // time step multiplier/divisor
double ls_percent=0.66; // threshold in per cent of max.iterations used to make a decision to change time step length
// soil parameters
// VGM parameters
double *vgm_ns,*vgm_s0s,*vgm_s1s,*vgm_as,*vgm_h0,*vgm_h1,*vgm_k,*vgm_power,*vgm_specific_storage; // VGM coefs per layer of soil
// thermal parameters
double *soil_lambda, *soil_Ct;
// diffusive parameters
double *soil_D[4];
int vgm_nlayers; // number of layers

///////////////////////////////////////
////////// basic solver class /////////
///////////////////////////////////////

class basic_solver;
class basic_solver {
public:
	// storage
	double *b_U; // solution (pressure head/concentration/temperature)
	double *sb_U; // to save solution
	double *MM,*RP; // multiplication results and right part for TFQMR
	// matrix parts
	double *pr_A=NULL,*pr_B=NULL,*pr_R=NULL;
	// soil parameters
	// VGM parameters
	double *pr_vgm_ns=NULL,*pr_vgm_s0s,*pr_vgm_s1s,*pr_vgm_as,*pr_vgm_k,*pr_vgm_power,*pr_vgm_specific_storage; // VGM coefs per cell
	double *pr_soil_lambda;
	double *pr_soil_Ct;
	double *pr_soil_D[4];

	// steps and auxiliary
	double tau; // current time step length 
	double T; // current time
	double L,dL; // domain depth and space variable step length
	int steady_state=0;
	// output
	FILE *out;
	// auxiliary
	// precalculate soil parameters for each cell in pr_*
	void soil_calc_coefs()
	{
	    if (pr_vgm_ns==NULL)
	    {
		pr_vgm_ns=new double[N+2];
		pr_vgm_s0s=new double[N+2];
		pr_vgm_s1s=new double[N+2];
		pr_vgm_as=new double[N+2];
		pr_vgm_k=new double[N+2];
		pr_vgm_power=new double[N+2];
		pr_vgm_specific_storage=new double[N+2];
		pr_soil_lambda=new double[N+2];
		pr_soil_Ct=new double[N+2];
		for (int i=0;i<4;i++)
			pr_soil_D[i]=new double[N+2];
		memset(pr_vgm_s0s,0,(N+2)*sizeof(double));
		// calculate
		for (int i=0;i<=N;i++)
		{
		double vgm_s0,vgm_s1,vgm_a,vgm_n,k,avr_power,vgm_ss,soil_l,soil_ct,sD[4];
	        avr_power=3.5;
		if (vgm_nlayers!=0)
		{
		    // i inside a layer
		    for (int j=0;j<vgm_nlayers;j++)
		    if (((i*dL)>=vgm_h0[j])&&((i*dL)<vgm_h1[j]))
		    {
			vgm_s0=vgm_s0s[j];
			vgm_s1=vgm_s1s[j];
			vgm_a=vgm_as[j];
			vgm_n=vgm_ns[j];
			vgm_ss=vgm_specific_storage[j];
			soil_l=soil_lambda[j];
			soil_ct=soil_Ct[j];
			for (int kk=0;kk<4;kk++) sD[kk]=soil_D[kk][j];
			k=vgm_k[j];
			avr_power=vgm_power[j];
			goto save;
		    }
		    double h1,h0,minh0=1e300,maxh1=0;
		    int i1,i0;
		    int mi1,mi0;
		    int first=1;
		    for (int j=0;j<vgm_nlayers;j++)
		    {
			if (vgm_h0[j]<minh0)
			{
			    minh0=vgm_h0[j];
			    mi0=j;
			}
			if (vgm_h1[j]>maxh1)
			{
			    maxh1=vgm_h1[j];
			    mi1=j;
			}
			for (int k=0;k<vgm_nlayers;k++)
			if (j!=k)
			if (((i*dL)>=vgm_h1[j])&&((i*dL)<vgm_h0[k]))
			{
			    if (first)
			    {
				h1=vgm_h1[j];
				i1=j;
				h0=vgm_h0[k];
				i0=k;
				first=0;
			    }
			    if (vgm_h1[j]>h1)
			    {
				h1=vgm_h1[j];
				i1=j;
			    }
			    if (vgm_h0[k]<h0)
			    {
				h0=vgm_h0[k];
				i0=k;
			    }
			}
		    }
		    // i between two layers
		    if (first==0)
		    {
			double _k=((i*dL)-h1)/(h0-h1);
			vgm_s0=_k*vgm_s0s[i0]+(1.0-_k)*vgm_s0s[i1];
			vgm_s1=_k*vgm_s1s[i0]+(1.0-_k)*vgm_s1s[i1];
			vgm_a=_k*vgm_as[i0]+(1.0-_k)*vgm_as[i1];
			vgm_n=_k*vgm_ns[i0]+(1.0-_k)*vgm_ns[i1];
			vgm_ss=_k*vgm_specific_storage[i0]+(1.0-_k)*vgm_specific_storage[i1];
			soil_l=_k*soil_lambda[i0]+(1.0-_k)*soil_lambda[i1];
			soil_ct=_k*soil_Ct[i0]+(1.0-_k)*soil_Ct[i1];
			for (int kk=0;kk<4;kk++) sD[kk]=_k*soil_D[kk][i0]+(1.0-_k)*soil_D[kk][i1];
			k=_k*vgm_k[i0]+(1.0-_k)*vgm_k[i1];
		    }
		    else
		    {
			// i below minimal
			if ((i*dL)<=minh0)
			{
			    vgm_s0=vgm_s0s[mi0];
			    vgm_s1=vgm_s1s[mi0];
			    vgm_a=vgm_as[mi0];
			    vgm_n=vgm_ns[mi0];
			    vgm_ss=vgm_specific_storage[mi0];
			    soil_l=soil_lambda[mi0];
			    soil_ct=soil_Ct[mi0];
			    for (int kk=0;kk<4;kk++) sD[kk]=soil_D[kk][mi0];
			    k=vgm_k[mi0];
			    avr_power=vgm_power[mi0];
			}
			else
			{
			    // i above maximal
			    if ((i*dL)>=maxh1)
			    {
				vgm_s0=vgm_s0s[mi1];
				vgm_s1=vgm_s1s[mi1];
				vgm_a=vgm_as[mi1];
				vgm_n=vgm_ns[mi1];
				vgm_ss=vgm_specific_storage[mi1];
				soil_l=soil_lambda[mi1];
				soil_ct=soil_Ct[mi1];
				for (int kk=0;kk<4;kk++) sD[kk]=soil_D[kk][mi1];
				k=vgm_k[mi1];
				avr_power=vgm_power[mi1];
			    }
			}
		    }
	    	}
		// save
save:
		pr_vgm_s0s[i]=vgm_s0;
		pr_vgm_s1s[i]=vgm_s1;
		pr_vgm_as[i]=vgm_a;
		pr_vgm_ns[i]=vgm_n;
		pr_vgm_specific_storage[i]=vgm_ss;
		pr_vgm_k[i]=k;
		pr_vgm_power[i]=avr_power;
		pr_soil_lambda[i]=soil_l;
		pr_soil_Ct[i]=soil_ct;
		for (int kk=0;kk<4;kk++) pr_soil_D[kk][i]=sD[kk];
		}
	    }
	}
	// linearly interpolate a value from a list of [T,V] pairs for currebt T 
	double linear_from_pair(double *TaT,double *TaF,int nTa,double __T=-1)
	{
		double ret=0.0;
		if (steady_state==1) return TaF[0];
		if (nTa != 0)
		{
			// calc linear conbination
			int i;
			double p = 0.0;
			double _T=T;
			if (__T!=-1) _T=__T;
			while (_T>TaT[nTa-1]) _T-=TaT[nTa-1];		
			for (i = 0;i<nTa;i++)
				if (TaT[i]>_T)
					break;
			if ((i == 0) || (i == nTa))
			{
				if (i == nTa) i--;
				p= TaF[i];
			}
			else
			{
				double k = (_T - TaT[i - 1]) / (TaT[i] - TaT[i - 1]);
				p= k*TaF[i] + (1 - k)*TaF[i - 1];
			}
			ret=p;
		}
		return ret;
	}
	// linear solver
	// linear system coefficients - A*U(i-1)+R*U(i)+B*U(i+1)=Rp
	virtual double A_(int i,double *pickard_prevU)=0;
	virtual double B(int i,double *pickard_prevU)=0;
	virtual double R(int i,double *pickard_prevU)=0;
	virtual double Rp(int i,double *U,double *pickard_prevU)=0;
	// matrix multiplication of vector UU
	virtual void mmult(double *UU)=0;
	void mmult_main(double *UU)
	{
#pragma omp for nowait
		for (int i = 1;i < N ;i++)
		{
			double rr=pr_R[i];
			MM[i]=UU[i-1]*pr_A[i]+UU[i+1]*pr_B[i]+UU[i]*rr;
			// normalization
			if (rr!=0.0)
				MM[i]/=rr;
		}
	}
	// precalculations for current time step
	virtual void __precalc(double *pickard_prevU)=0;
	void __precalc_main(double *pickard_prevU)
	{
#pragma omp for nowait
		for (int i = 0;i <= N ;i++)
		{
			pr_A[i]=A_(i,pickard_prevU);
			pr_B[i]=B(i,pickard_prevU);
			pr_R[i]=R(i,pickard_prevU);
			// right part
			RP[i]=Rp(i,b_U,pickard_prevU);
			// normalization
			if (pr_R[i]!=0.0)
				RP[i]/=pr_R[i];
		}
	}
	// linear solver (TFQMR) (returns 1 if tau was increased at the end)
	int calc_step(double *pickard_prevU=NULL,int fix_tau=0)
	{
		if (pickard_prevU==NULL) pickard_prevU=b_U;
		// first step
		if (T==tau)
		{
		    for (int i =0;i <= N+1 ;i++)
		    {
			MM[i]=0.0;
			RP[i]=0.0;
		    }
		}
		//////////
		double *w=new double[(N+2)];
		double *y[2];
		y[0]=new double[(N+2)];
		y[1]=new double[(N+2)];
		double *rr=new double[(N+2)];
		double *v=new double[(N+2)];
		double *d=new double[(N+2)];
		double *x=new double[(N+2)];
		double theta=0,nu=0,tt=0,ro=0,ro1=0,c=0,b=0;
		double sgm,aa,rv;
		int n,m;
start:
		theta=nu=tt=ro=ro1=c=b=0.0;
		memcpy(x,b_U,(N+2)*sizeof(double));
#pragma omp parallel
	    {
		////////// precalculations
		__precalc(pickard_prevU);
		mmult(x);
#pragma omp for reduction(+:tt)
		for (int i=0;i<(N+2);i++)
		{
			// w_1=y_1=r_0=b-Ax_0
			// rr0: ro=rr_0*r_0!=0 - rr0=r0
			w[i]=y[0][i]=RP[i]-MM[i];
			// d=0
			d[i]=0;
			// tau=||r0||
			tt+=w[i]*w[i];
		}
#pragma omp single
		tt=sqrt(tt);
#pragma omp barrier
		// random rr0, ro0=rr_0*r_0
#pragma omp for reduction(+:ro)
		for (int i=0;i<(N+2);i++)
		{
			rr[i]=tt*((rand() % 10000) / (10000.0 - 1.0));
			ro+=rr[i]*w[i];
		}
		// v=Ay_1
		mmult(y[0]);
#pragma omp for
		for (int i=0;i<(N+2);i++)
			v[i]=MM[i];
#pragma omp single
		n=1;
#pragma omp barrier
loop1:
		{
			int br=0;
			// sigma_n-1 - rr_0*v_n-1
#pragma omp single
			sgm=0;
#pragma omp barrier
#pragma omp for reduction(+:sgm)
			for (int i=0;i<(N+2);i++)
				sgm+=rr[i]*v[i];
			// a_n-1=po_n-1/sigma_n-1
#pragma omp single
			aa=ro/sgm;
#pragma omp barrier
			// y_2n=y_2n-1 - a_n-1 * v_n-1
#pragma omp for
			for (int i=0;i<(N+2);i++)
				y[1][i]=y[0][i]-aa*v[i];
#pragma omp single
			m=2*n-1;
#pragma omp barrier
loop2:
			{
				double ot=theta,onu=nu;
				// w_m+1=w_m-a_n-1Ay_m
				mmult(y[m-(2*n-1)]);
#pragma omp for
				for (int i=0;i<(N+2);i++)
					w[i]=w[i]-aa*MM[i];
				// theta_m=||w_m+1||/tau_m-1; c_m=1/sqrt(1+theta_m^2)
#pragma omp single
				theta=0;
#pragma omp barrier
#pragma omp for reduction(+:theta)
				for (int i=0;i<(N+2);i++)
					theta+=w[i]*w[i];
#pragma omp single
			    {
				theta=sqrt(theta)/tt;
				c=1.0/sqrt(1.0+theta*theta);
				// tau_m=tau_m-1 * theta_m *c_m; nu_m=c_m^2 *a_n-1
				tt=tt*theta*c;
				nu=c*c*aa;
				rv=0.0;
			    }
#pragma omp barrier
				// d_m = y_m+(theta_m-1^2 nu_m-1 / a_n-1)*d_m-1
				// x_m=x_m-1 + nu_m *d_m
#pragma omp for
				for (int i=0;i<(N+2);i++)
				{
					d[i]=y[m-(2*n-1)][i]+d[i]*(ot*ot*onu/aa);
					x[i]=x[i]+nu*d[i];
				}
				mmult(x);
#pragma omp for reduction(+:rv)
				for (int i=0;i<(N+2);i++)
					rv+=(RP[i]-MM[i])*(RP[i]-MM[i]);
#pragma omp single
				rv=sqrt(rv)/((N+2));
#pragma omp barrier
				if (rv<ls_eps)
				{
				    br=1;
				    goto eloop2;
				}
			}
#pragma omp single
			m++;
#pragma omp barrier
			if (m<=2*n)
			    goto loop2;
eloop2:
			if (br==1)
				goto eloop1;
			// ro_n=rr0*w_2n+1, beta_n=ro_n/ro_n-1
#pragma omp single
			ro1=0;
#pragma omp barrier
#pragma omp for reduction(+:ro1)
			for (int i=0;i<(N+2);i++)
				ro1+=rr[i]*w[i];
#pragma omp single
		    {
			b=ro1/ro;
			ro=ro1;
		    }
#pragma omp barrier
			// y_2n+1 = w_2n+1+beta_n*y_2n
#pragma omp for
			for (int i=0;i<(N+2);i++)
				y[0][i]=w[i]+b*y[1][i];
			// v_n=Ay_2n+1+b*(Ay_2n+b*v_n-1)
			mmult(y[1]);
#pragma omp for
			for (int i=0;i<(N+2);i++)
				v[i]=b*(MM[i]+b*v[i]);
			mmult(y[0]);
#pragma omp for
			for (int i=0;i<(N+2);i++)
				v[i]=MM[i]+v[i];
		}
#pragma omp single
		n++;
#pragma omp barrier
		if (n<ls_max_iter)
		    goto loop1;
eloop1:;
	    }
		// change tau and recalc if needed
		if (fix_tau==0)
		if (n==ls_max_iter)
		{
		     if (tau<ls_min_tau)
		     {
			 x[0]=NAN;
			 printf("minimal tau value reached\n");
			 exit(2);
		     }
		     else
		     {
		    	T-=tau;
		    	tau/=ls_mult;
		    	T+=tau;
			// debug output
			if (debug_level==2) printf("%p r n %d r %g T %g tau %g sol %g %g %g %g %g %g\n",this,n,rv,T,tau,x[0],x[1],x[2],b_U[0],b_U[1],b_U[2]);
		        goto start;
		     }
		}
		// save solution and free memory
		memcpy(b_U,x,(N+2)*sizeof(double));
		delete [] w;
		delete [] y[0];
		delete [] y[1];
		delete [] rr;
		delete [] v;
		delete [] d;
		delete [] x;
		// increase tau if needed and increase T
		int ret=0;
		if (fix_tau==0)
		if (n<ls_max_iter*ls_percent)
			if (tau<max_tau)
			{
				 tau*=ls_mult;
				 ret=1;
			}
		T+=tau;
		// debug output
		if (debug_level==2) printf("%p e n %d rv %g T %g tau %g\n",this,n,rv,T,tau);
		return ret;
	}
	// perform Pickard iteration 
	void pickard_calc_step(basic_solver **other=NULL,int nother=0)
	{
	    double *prev_U=new double[N+2];
	    double *init_U=new double[N+2];
	    double **prev_other_U=NULL,**init_other_U=NULL;
	    double diff=0;
	    int iter=0;
	    int mul=0;
	    int *other_muls=NULL;
	    memcpy(init_U,b_U,sizeof(double)*(N+2));
restart:
	    memcpy(b_U,init_U,sizeof(double)*(N+2));
	    if (init_other_U)
	    for (int i=0;i<nother;i++)
		memcpy(other[i]->b_U,init_other_U[i],sizeof(double)*(N+2));
	    iter=0;
	    mul=calc_step();
	    if (mul)
	    {	
		T-=tau;
		tau/=ls_mult;
		T+=tau;
	    }
	    if (nother) // run other solvers
	    {
		if (prev_other_U==NULL) { prev_other_U=new double *[nother]; memset(prev_other_U,0,nother*sizeof(void *));}
		if (init_other_U==NULL) { init_other_U=new double *[nother]; memset(init_other_U,0,nother*sizeof(void *));}
		if (other_muls==NULL) other_muls=new int[nother];
		for (int i=0;i<nother;i++)
		{
			if (prev_other_U[i]==NULL) prev_other_U[i]=new double [N+2];
			if (init_other_U[i]==NULL) init_other_U[i]=new double [N+2];
			memcpy(init_other_U[i],other[i]->b_U,sizeof(double)*(N+2));
			other_muls[i]=other[i]->calc_step();
			if (other_muls[i])
			{	
				other[i]->T-=other[i]->tau;
				other[i]->tau/=ls_mult;
				other[i]->T+=other[i]->tau;
		        }
		}
		// find minimal tau
		double min_tau=tau;
		for (int i=0;i<nother;i++)
			if (other[i]->tau<min_tau)
				min_tau=other[i]->tau;
		// solve all once more with fixed minimal tau
		memcpy(b_U,init_U,sizeof(double)*(N+2));
		for (int i=0;i<nother;i++)
			memcpy(other[i]->b_U,init_other_U[i],sizeof(double)*(N+2));
		T-=2.0*tau;
		tau=min_tau;
		T+=tau;
		mul=calc_step(NULL,1);
		for (int i=0;i<nother;i++)
		{
			other[i]->T-=2.0*other[i]->tau;
			other[i]->tau=min_tau;
			other[i]->T+=other[i]->tau;
			other[i]->calc_step(NULL,1);
		}
	    }
	    // debug
	    if (debug_level==1)
	    {
		 printf("%p pickard init tau %g T %g other %d ",this,tau,T-tau,nother);
		for (int i=0;i<nother;i++) printf(" %d tau %g T %g ",i,other[i]->tau,other[i]->T-other[i]->tau);
		printf("\n");
	    }
	    do
	    {
		// save solution on previous iteration and restore initial values
		memcpy(prev_U,b_U,sizeof(double)*(N+2));
		memcpy(b_U,init_U,sizeof(double)*(N+2));
		for (int i=0;i<nother;i++)
		{
			memcpy(prev_other_U[i],other[i]->b_U,sizeof(double)*(N+2));
			memcpy(other[i]->b_U,init_other_U[i],sizeof(double)*(N+2));
		}
		// solve on next iteration
		T-=tau;
		calc_step(prev_U,1);
		for (int i=0;i<nother;i++)
		{
		    other[i]->T-=other[i]->tau;
		    other[i]->calc_step(prev_other_U[i],1);
		}
		// calculate difference
		diff=0;
		for (int i=0;i<N+1;i++)
		    diff+=(b_U[i]-prev_U[i])*(b_U[i]-prev_U[i]);
		for (int j=0;j<nother;j++)
		    for (int i=0;i<N+1;i++)
			diff+=(other[j]->b_U[i]-prev_other_U[j][i])*(other[j]->b_U[i]-prev_other_U[j][i]);
		diff=diff/(N*(nother?nother:1));
		// debug
		if (debug_level==1) printf("%p pickard iter %d diff %g T %g S %g %g %g - %g %g %g\n",this,iter,diff,T-tau,b_U[0],b_U[1],b_U[2],b_U[N-2],b_U[N-1],b_U[N]);
		if ((iter++)> ls_max_iter) 
			break;
	    }
	    while (diff>ls_eps);
	    if (!isfinite(diff))
	    {
		ls_max_iter/=2.0;
		goto restart;
	    }
	    if (iter>=ls_max_iter) // decrease step
	    {
		 if (tau<ls_min_tau)
		 {
			 b_U[0]=NAN;
			 printf("minimal tau value reached\n");
			 exit(2);
		 }
		 else
		 {
		    T-=tau;
		    tau/=ls_mult;
		    T+=tau;
		    if (debug_level==1) printf("%p pickard restart niter %d diff %g tau %g T %g\n",this,iter,diff,tau,T-tau);
		    goto restart;
		}
	    }
	    // increase tau if needed
	    if (mul)
	    {
		T-=tau;
		tau*=ls_mult;
		T+=tau;
	    }
	    if (nother)
		for (int i=0;i<nother;i++)
		    if (other_muls[i])
		    {
			other[i]->T-=other[i]->tau;
			other[i]->tau*=ls_mult;
			other[i]->T+=other[i]->tau;
		    }
	    // debug
	    if (debug_level==1) printf("%p pickard niter %d diff %g tau %g T %g\n",this,iter,diff,tau,T-tau);
	    // clean up
	    delete [] prev_U;
	    delete [] init_U;
	    if (nother)
	    {
		for (int i=0;i<nother;i++)
		{
		    delete [] prev_other_U[i];
		    delete [] init_other_U[i];
		}
		delete [] prev_other_U;
		delete [] init_other_U;
		delete [] other_muls;
	    }
	}
	// solve up to the time point eT
	void solve_up_to_T(double eT,double _tau,basic_solver **other=NULL,int nother=0)
	{
	    while ((T-tau)<(eT-ls_min_tau))
	    {
		if (T+tau>eT)
		{
			    T-=tau;
			    tau=eT-T;
			    T+=tau;
		}
		pickard_calc_step(other,nother);
	    }
	}
	// constructor
	basic_solver(double _L)
	{
		b_U=new double[(N+2)];
		sb_U=new double[(N+2)];
		MM=new double[(N+2)];
		RP=new double[(N+2)];

		L = _L;
		dL = L/N;
		tau = max_tau;
		T=tau;
	}
	// desctructor
	~basic_solver()
	{
		delete [] b_U;
		delete [] sb_U;
		delete [] MM;
		delete [] RP;
		if (pr_A) delete [] pr_A;
		if (pr_B) delete [] pr_B;
		if (pr_R) delete [] pr_R;
	}
	// output
	virtual void output_solution()=0;
};
//////////////////////////////////////
// Solver for Richards equation
//////////////////////////////////////
class H_solver: public basic_solver {
public:
	double *b_V; // water movement velocity
	double *b_Vd; // d velocity / dz
	// auxiliary
	double Dg; // auxiliary
	double H0_; // initial H parameter
	FILE *v_out;
	// inputs
	// EV per day
	double *EV_T; // times
	double *EV_F; // values
	int nEV_T; // size of EV_T and EV_F
	double pr_ET; // current evapotranspiration
	// precipitation
	double *perc_T, *perc_A; // times and values
	int nperc; // size
	double pr_perc; // current precipitation
	// grownwater level
	double *flT, *flV; // times and values
	int nfl; // size of flT,flV
	double pr_fl; // precalculated growndwater level

	// values precalculated on the start of time step
	double *pr_dwdh=NULL; // d theta / dh
	double *pr_w=NULL; // wetness
	double *pr_K=NULL; // hydraulic conductivity

	/// calc wetness on the base of van Genuchten model ////////
	double wetness(int i,double P)
	{
		P *= HtoPmult;
		if (P <= 0.0)
			return pr_vgm_s1s[i];
		return pr_vgm_s0s[i]+((pr_vgm_s1s[i] - pr_vgm_s0s[i]) / pow(1 + pow(pr_vgm_as[i]*P*10.19716, pr_vgm_ns[i]), (1 - 1 / pr_vgm_ns[i])));
	}
	// precalulculate wetness in all cells in pr_w
	void pr_wetness()
	{
	    // alloc
#pragma omp critical
	    if (pr_w==NULL)
		pr_w=new double[N+2];
#pragma omp barrier
#pragma omp for
	    for (int i=0;i<=N;i++)
		pr_w[i]=wetness(i,-b_U[i]);
#pragma omp barrier
	}
	// precalculate d theta / dh in all cells in pr_dwdh
	double dw_dh(double P,int i)
	{
		double Ch;
		if (P <= 0.0)
			Ch = 0.0;
		else
		{
			P*=10.19716;
			P*=HtoPmult;
			double aPn=pow(pr_vgm_as[i]*P, pr_vgm_ns[i]);
			Ch=-(1.0/(P/(HtoPmult*10.197196)))*(((1.0/pr_vgm_ns[i])-1)*pr_vgm_ns[i]*aPn*pow(aPn+1.0,(1.0/pr_vgm_ns[i])-2.0)*(pr_vgm_s1s[i]-pr_vgm_s0s[i]));
		}
		return Ch;
	}
	void pr_dw_dh()
	{
	    // alloc
#pragma omp critical
	    if (pr_dwdh==NULL)
		pr_dwdh=new double[N+2];
#pragma omp barrier
#pragma omp for
	    for (int i=0;i<=N;i++)
	    {
		double P;
		double Ch = 0.0;
		int idx=0;
		P = -b_U[idx=i];
		// calculate
		Ch = dw_dh(P,i);
		if (!isfinite(Ch)) Ch=0;
		pr_dwdh[idx]=Ch;
	    }
#pragma omp barrier
	}
	double Kr(double w,int i)
	{
		// Mualem's model
		double kr = pow(((w - pr_vgm_s0s[i]) / (pr_vgm_s1s[i] - pr_vgm_s0s[i])), -pr_vgm_power[i])
		    *pow(1.0 - pow(1.0 - pow(((w - pr_vgm_s0s[i]) / (pr_vgm_s1s[i] - pr_vgm_s0s[i])), 
			    (1.0 / (1.0 - 1.0 / pr_vgm_ns[i]))), (1.0 - 1.0 / pr_vgm_ns[i])), 2.0);
		return kr;
	}
	// precalculate hydraulic conductivity in all cells in pr_K
	void pr_KK()
	{
	    // alloc
#pragma omp critical
	    if (pr_K==NULL)
		pr_K=new double[N+2];
#pragma omp barrier
#pragma omp for
	    for (int i=0;i<=N;i++)
		pr_K[i]=Kr(pr_w[i],i)*pr_vgm_k[i];
#pragma omp barrier
	}
	// sets current ET, precipitation, root system depth, groundwater level
	void precalc_values()
	{
		double v=0.0,k;
		int i=0;
		// precipitation
		pr_perc=linear_from_pair(perc_T,perc_A,nperc);
		// Evapotranspiration
		pr_ET=linear_from_pair(EV_T,EV_F,nEV_T);
		// grownwater level
		pr_fl=linear_from_pair(flT,flV,nfl);

		if (pr_fl>0) pr_perc=pr_ET=0; // no flow if soil is submerged under water
	}
	// testing
	double testing_f(double z,double t) // solution for testing mode
	{
		return t*t*(z*z-1)/(30.0*86400*30.0*86400);
	}
	double testing_dfdt(double z,double t)
	{
		return  2*t*(z*z-1)/(30.0*86400*30.0*86400);
	}
	double testing_dfdz(double z,double t)
	{
		return  2*z*t*t/(30.0*86400*30.0*86400);
	}
	double testing_d2fdz2(double z,double t)
	{
		return  2*t*t/(30.0*86400*30.0*86400);
	}
	double testing_dkdh(double z,double t)
	{
		double S=(wetness(0,-testing_f(z,t))-pr_vgm_s0s[0])/(pr_vgm_s1s[0]-pr_vgm_s0s[0]);
		double dSdh=dw_dh(-testing_f(z,t),0)/(pr_vgm_s1s[0]-pr_vgm_s0s[0]);
		
		if (S==1.0) return 0.0;
		double ret=-pr_vgm_power[0]*pow(S,-pr_vgm_power[0]-1)*
			    pow(1.0 - pow(1.0 - pow(S, (1.0 / (1.0 - 1.0 / pr_vgm_ns[0]))), (1.0 - 1.0 / pr_vgm_ns[0])), 2.0);
		ret+=2.0*pow(S,-pr_vgm_power[0])*(1.0 - pow(1.0 - pow(S, (1.0 / (1.0 - 1.0 / pr_vgm_ns[0]))), (1.0 - 1.0 / pr_vgm_ns[0])))*
					pow(1.0 - pow(S, (1.0 / (1.0 - 1.0 / pr_vgm_ns[0]))), (1.0 - 1.0 / pr_vgm_ns[0])-1)*
					pow(S, (1.0 / (1.0 - 1.0 / pr_vgm_ns[0]))-1);
		ret*=pr_vgm_k[0]*dSdh;
		return ret;
	}
	double testing_F(double z,double t)
	{
		double F;
		F=(dw_dh(-testing_f(z,t),0)+wetness(0,-testing_f(z,t))*pr_vgm_specific_storage[0]/pr_vgm_s1s[0])*testing_dfdt(z,t);
		F-=testing_dkdh(z,t)*testing_dfdz(z,t)*(testing_dfdz(z,t)-1.0)+pr_vgm_k[0]*Kr(wetness(0,-testing_f(z,t)),0)*testing_d2fdz2(z,t);
		return F;
	}
	// upper boundary condition (DbU=Uc)
	double Uc()
	{
	    if (testing_mode==1)
		return dL*testing_dfdz(0,T);
	    // condition for flux
	    double kk=(pr_K[1]+pr_K[0])/2;
	    double ret=dL*(((pr_ET-pr_perc)/kk)+1);
	    return ret;
	}
	// coefficients of three-diagonal linear equations system
	double Km1(int ii,double *pickard_prevU)
	{
		double km1;
		if (ii!=0) 
		    km1=pr_K[ii-1];
		else
		    km1=pr_vgm_k[0]*Kr(wetness(0,-(pickard_prevU[1]-2.0*Uc())),0);
		return km1;
	}
	double A_(int i,double *pickard_prevU)
	{
		return Dg*(0.25*pr_K[i+1]-0.25*Km1(i,pickard_prevU)-pr_K[i]);
	}
	double B(int i,double *pickard_prevU)
	{
		return Dg*(-0.25*pr_K[i+1]+0.25*Km1(i,pickard_prevU)-pr_K[i]);
	}
	double R(int i,double *pickard_prevU)
	{
		return 2.0*pr_K[i]*Dg+((steady_state==0)?((pr_dwdh[i]+(pr_w[i] / pr_vgm_s1s[i])*pr_vgm_specific_storage[i])/tau):0);
	}
	// right part 
	double Rp(int ii,double *b_Uold,double *pickard_prevU)
	{
		double ret=0.0;
		ret = ((steady_state==0)?(((pr_dwdh[ii]+(pr_w[ii] / pr_vgm_s1s[ii])*pr_vgm_specific_storage[ii])/tau)*b_Uold[ii]):0);
		ret-=(pr_K[ii+1]-Km1(ii,pickard_prevU))/(2.0*dL);
		if (testing_mode==1) ret+=testing_F(ii*dL,T);
		return ret;
	}
	// linear solver
	// matrix multiplication of vector UU
	void mmult(double *UU)
	{
		mmult_main(UU);
#pragma omp single
		{
			// bottom boundary condition - h level
			MM[N]=UU[N];
			// upper boundary condition - flux
			MM[0]=UU[0]+((pr_A[0]+pr_B[0])/pr_R[0])*UU[1]; // second order condition through "ghost point"
		}
#pragma omp barrier
	}
	// precalculations for current time step
	void __precalc(double *pickard_prevU)
	{
#pragma omp single
		{ // precalculate non-linear coefficients using the solution of previous iteration
		    memcpy(sb_U,b_U,sizeof(double)*(N+2));
		    memcpy(b_U,pickard_prevU,sizeof(double)*(N+2));
		}
#pragma omp barrier
		pr_wetness();
		pr_dw_dh();
		pr_KK();
		// precalc matrix coefficients
#pragma omp single
	    {
		precalc_values();
		if (pr_A==NULL) pr_A=new double[N+2];
		if (pr_B==NULL) pr_B=new double[N+2];
		if (pr_R==NULL) pr_R=new double[N+2];
		memcpy(b_U,sb_U,sizeof(double)*(N+2));
	    }
#pragma omp barrier
	    // precalc linear system values
		__precalc_main(pickard_prevU);
#pragma omp barrier
#pragma omp single
		{
			// bottom boundary condition - groundwater level through given h
			RP[N]=L+pr_fl;
			// upper boundary condition - flux
			RP[0]=(Rp(0,b_U,pickard_prevU)+pr_A[0]*2.0*Uc())/pr_R[0];  // second order condition through "ghost point"
			if (testing_mode==1)
			    RP[N]=testing_f(N*dL,T);
		}
#pragma omp barrier
	}
	void calculate_velocity()
	{
		for (int i=1;i<N;i++)
		    b_V[i]=pr_K[i]*(((b_U[i+1]-b_U[i-1])/(2.0*dL))-1.0);
		b_V[0]=pr_K[0]*(((b_U[1]-b_U[0])/dL)-1);
		b_V[N]=pr_K[N]*(((b_U[N]-b_U[N-1])/dL)-1);
		for (int i=1;i<N;i++)
		    b_Vd[i]=(b_V[i+1]-b_V[i-1])/(2.0*dL);
		b_Vd[0]=(b_V[1]-b_V[0])/dL;
		b_Vd[N]=(b_V[N]-b_V[N-1])/dL;

		if (testing_mode==2) // no advection
		    for (int i=0;i<=N;i++)
			b_V[i]=b_Vd[i]=0.0;
	}
	// constructor
	H_solver(double _L,double H0) : basic_solver(_L)
	{
		b_V=new double [N+2];
		b_Vd=new double [N+2];
		Dg=pow(dL,-2);
		H0_=H0;

		// initial conditions
		for (int i = 0;i < N + 1;i++)
		{
			b_U[i] = H0+i*dL;
			if (testing_mode==1) b_U[i]=testing_f(i*dL,0);
		}
		out=fopen("out_H.txt","wt");
		v_out=fopen("out_V.txt","wt");
	}
	// desctructor
	~H_solver()
	{
		delete b_V;
		delete b_Vd;
		if (pr_K) delete [] pr_K;
		if (pr_w) delete [] pr_w;
		if (pr_dwdh) delete [] pr_dwdh;
		if (vgm_ns) delete [] vgm_ns;
		if (vgm_s0s) delete [] vgm_s0s;
		if (vgm_s1s) delete [] vgm_s1s;
		if (vgm_as) delete [] vgm_as;
		if (vgm_h0) delete [] vgm_h0;
		if (vgm_h1) delete [] vgm_h1;
		if (vgm_k) delete [] vgm_k;
		if (vgm_power) delete [] vgm_power;
		if (vgm_specific_storage) delete [] vgm_specific_storage;
		if (pr_vgm_ns) delete [] pr_vgm_ns;
		if (pr_vgm_s0s) delete [] pr_vgm_s0s;
		if (pr_vgm_s1s) delete [] pr_vgm_s1s;
		if (pr_vgm_as) delete [] pr_vgm_as;
		if (pr_vgm_k) delete [] pr_vgm_k;
		if (pr_vgm_power) delete [] pr_vgm_power;
		if (pr_vgm_specific_storage) delete [] pr_vgm_specific_storage;
		if (EV_T) delete [] EV_T;
		if (EV_F) delete [] EV_F;
		if (perc_T) delete [] perc_T;
		if (perc_A) delete [] perc_A;
		if (flT) delete [] flT;
		if (flV) delete [] flV;
		fclose(out);
		fclose(v_out);
	}
	// output
	void output_solution()
	{
		if (pr_K) calculate_velocity();
		fprintf(out,"t(days) %g tau(seconds) %g ET %g Prec %g Fl %g - H: ",(T-tau)/86400.0,tau,pr_ET,pr_perc,pr_fl);
		if (T==tau)
		{
		    for (int i=0;i<N+1;i+=zoutstep)
			fprintf(out,"%g ",i*dL);
		    fprintf(out,"\n");
		    fprintf(out,"t(days) %g tau(seconds) %g ET %g Prec %g Fl %g - H: ",(T-tau)/86400.0,tau,pr_ET,pr_perc,pr_fl);
		}
		for (int i=0;i<N+1;i+=zoutstep)
			fprintf(out,"%g ",b_U[i]);
		// calculate water table value (saturated zone starting from the bottom)
		double wt;
		for (wt=N;wt>0;wt--)
			if (b_U[(int)wt]<0)
				break;
		wt*=dL;
		fprintf(out,"water_table %g ",wt);
		// velocities
		fprintf(v_out,"t(days) %g V: ",(T-tau)/86400.0);
		if (T==tau)
		{
		    for (int i=0;i<N+1;i+=zoutstep)
			fprintf(v_out,"%g ",i*dL);
		    fprintf(v_out,"\n");
		    fprintf(v_out,"t(days) %g V: ",(T-tau)/86400.0);
		}
		for (int i=0;i<N+1;i+=zoutstep)
			fprintf(v_out,"%g ",b_V[i]);
		fprintf(v_out,"\n");
		if (testing_mode==1)
		{
		    double diff=0,sq1=0,sq2=0;
		    for (int i=0;i<N+1;i++)
		    {
			diff+=(b_U[i]-testing_f(i*dL,T-tau))*(b_U[i]-testing_f(i*dL,T-tau));
			sq1+=b_U[i]*b_U[i];
			sq2+=testing_f(i*dL,T-tau)*testing_f(i*dL,T-tau);
		    }
		    diff=sqrt(diff)/N;
		    sq1=sqrt(sq1);
		    sq2=sqrt(sq2);
		    printf ("testing diff %g rel %g\n",diff,fabs(sq1-sq2)/sq2);
		}
		fprintf(out,"\n");
	}
};
///////////////////////////////////////////////////////
// Solver for temperature transport  equation ////////
///////////////////////////////////////////////////////
class T_solver: public basic_solver {
public:
	// inputs
	// atmospheric temperature
	double *TaT, *TaF; // times and values
	int nTa; // size
	double pr_Ta; // precalculated value
	H_solver *Hs; // solver for H equation

	// sets current Ta
	void precalc_values()
	{
		pr_Ta=linear_from_pair(TaT,TaF,nTa);
	}
	// testing
	double testing_f(double z,double t) // solution for testing mode
	{
		return t*t*z*z/(30.0*86400*30.0*86400);
	}
	double testing_dfdt(double z,double t)
	{
		return  2*t*z*z/(30.0*86400*30.0*86400);
	}
	double testing_dfdz(double z,double t)
	{
		return  2*z*t*t/(30.0*86400*30.0*86400);
	}
	double testing_d2fdz2(double z,double t)
	{
		return  2*t*t/(30.0*86400*30.0*86400);
	}
	double testing_F(double z,double t)
	{
		return pr_soil_Ct[0]*testing_dfdt(z,t)-pr_soil_lambda[0]*testing_d2fdz2(z,t)+pr_soil_Ct[0]*Hs->pr_vgm_k[0]*
		    Hs->Kr(Hs->wetness(0,-Hs->testing_f(z,t)),0)*(Hs->testing_dfdz(z,t)-1)*testing_dfdz(z,t)
		    +pr_soil_Ct[0]*(Hs->pr_vgm_k[0]*Hs->Kr(Hs->wetness(0,-Hs->testing_f(z,t)),0)*Hs->testing_d2fdz2(z,t)+
			    Hs->testing_dkdh(z,t)*Hs->testing_dfdz(z,t)*(Hs->testing_dfdz(z,t)-1))*testing_f(z,t);
	}
	// coefficients of three-diagonal linear equations system
	double A_(int i,double *pickard_prevU)
	{
		return -(pr_soil_lambda[i]/(dL*dL))-pr_soil_Ct[i]*Hs->b_V[i]/(2.0*dL);
	}
	double B(int i,double *pickard_prevU)
	{
		return -pr_soil_lambda[i]/(dL*dL)+pr_soil_Ct[i]*Hs->b_V[i]/(2.0*dL);
	}
	double R(int i,double *pickard_prevU)
	{
		return 2.0*pr_soil_lambda[i]/(dL*dL)+(pr_soil_Ct[i]*Hs->b_Vd[i])+((steady_state==0)?(pr_soil_Ct[i]/tau):0.0);
	}
	// right part 
	double Rp(int ii,double *b_Uold,double *pickard_prevU)
	{
		double ret=((steady_state==0)?(pr_soil_Ct[ii]/tau)*b_Uold[ii]:0.0);
		if (testing_mode==1)
		    ret+=testing_F(ii*dL,T);
		return ret;
	}
	// linear solver
	// matrix multiplication of vector UU
	void mmult(double *UU)
	{
		mmult_main(UU);
#pragma omp single
		{
			// bottom boundary condition - transparent flow
			MM[N]=UU[N]+((pr_A[N]+pr_B[N])/(pr_R[N]))*UU[N-1];  // second order condition through "ghost point"
			// upper boundary condition
			MM[0]=UU[0];
		}
	}
	// precalculations for current time step
	void __precalc(double *pickard_prevU)
	{
#pragma omp single
		{ // precalculate non-linear coefficients using the solution of previous iteration
		    memcpy(sb_U,b_U,sizeof(double)*(N+2));
		    memcpy(b_U,pickard_prevU,sizeof(double)*(N+2));
		    precalc_values();
		    if (pr_A==NULL) pr_A=new double[N+2];
		    if (pr_B==NULL) pr_B=new double[N+2];
		    if (pr_R==NULL) pr_R=new double[N+2];
		    memcpy(b_U,sb_U,sizeof(double)*(N+2));
		}
#pragma omp barrier
	    // precalc linear system values
		__precalc_main(pickard_prevU);
#pragma omp barrier
#pragma omp single
		{
			// bottom boundary condition - transparent flow
			RP[N]=Rp(N,b_U,pickard_prevU)/pr_R[N];  // second order condition through "ghost point"
			// upper boundary condition
			RP[0]=pr_Ta;
			if (testing_mode==1)
			{
				RP[0]=testing_f(0,T);
				RP[N]=(Rp(N,b_U,pickard_prevU)-pr_B[N]*2.0*dL*testing_dfdz(N*dL,T))/pr_R[N]; 
			}
		}
	}
	// constructor
	T_solver(double _L,H_solver *_Hs) : basic_solver(_L), Hs(_Hs)
	{
		out=fopen("out_T.txt","wt");
	}
	void initialize() // must be called after setting Ta
	{
		precalc_values();
		// initial conditions
		for (int i = 0;i < N + 1;i++)
		{
			b_U[i] = pr_Ta;
			if (testing_mode==1)
			    b_U[i]=testing_f(i*dL,0);
		}
	}
	// desctructor
	~T_solver()
	{
		if (TaT) delete [] TaT;
		if (TaF) delete [] TaF;
		fclose(out);
	}
	// output
	void output_solution()
	{
		fprintf(out,"t(days) %g tau(seconds) %g T: ",(T-tau)/86400.0,tau);
		if (T==tau)
		{
		    for (int i=0;i<N+1;i+=zoutstep)
			fprintf(out,"%g ",i*dL);
		    fprintf(out,"\n");
		    fprintf(out,"t(days) %g tau(seconds) %g T: ",(T-tau)/86400.0,tau);
		}
		for (int i=0;i<N+1;i+=zoutstep)
		    fprintf(out,"%g ",b_U[i]);
		if (testing_mode==1)
		{
		    double diff=0,sq1=0,sq2=0;
		    for (int i=0;i<N+1;i++)
		    {
			diff+=(b_U[i]-testing_f(i*dL,T-tau))*(b_U[i]-testing_f(i*dL,T-tau));
			sq1+=b_U[i]*b_U[i];
			sq2+=testing_f(i*dL,T-tau)*testing_f(i*dL,T-tau);
		    }
		    diff=sqrt(diff)/N;
		    sq1=sqrt(sq1);
		    sq2=sqrt(sq2);
		    printf ("testing diff %g rel %g\n",diff,fabs(sq1-sq2)/sq2);
		}
		fprintf(out,"\n");
	}
};
///////////////////////////////////////////////////////
// Solver for SOC compounds transport /////////////////
///////////////////////////////////////////////////////
class C_solver: public basic_solver {
public:
	// diffusion equation parameters
	double q; // fractional derivative order
	double fanox;
	double inputs_exponent;
	double inputs_divisor;
	double _gamma;
	int corrector_type=0; // type for the corrector for fractal dimensions
	double pr_corr;
	int compound_id=0; // compound id (0,1,2,3)
	double C0;
	// auxiliary
	double G2q,G1q;
	double prev_out=0;
	double balance,sum0;
	FILE *rho_out,*soc_out,*co2_out,*ch4_out;
	// RothC model parameters
	double alpha=0.1947; // for 24.5% clay content
	double beta=0.2371;
	double T0;
	double k[4]={0.8333,0.025,0.055,0.0017}; // monthly decomposition rate
	// oxydation parameters
	double ox,ox_tau;
	// storage for previous solutions needed to calculate fractional derviative's values
	std::vector<double *> prevC;
	std::vector<double> Ts; // moments of time
	// storage for previously calculated fluxes
	std::vector<double> Is,R0s,R1s,Os; // input, co2 respiration, ch4 respiration, outflow
	std::vector<double> fTs; // moments of time

	H_solver *Hs; // solver for H equation
	T_solver *Tsolv; // solver for T equation
	basic_solver *Css[4]; // solvers for other 3 C equations
	// inputs
	// SOC inputs
	double *SoCT, *SoCF; // times and values
	int nSoC; // size
	double pr_SoC; // precalculated value
	// soil cover coefficient
	double *KcT, *KcF; // times and values
	int nKc; // size
	double pr_Kc; // precalculated value

	// save to storage and clean up (finite memory)
	void save_to_storage()
	{
	    Ts.push_back(T-tau);
	    double *c=new double [N+2];
	    memcpy(c,b_U,sizeof(double)*(N+2));
	    prevC.push_back(c);
	    // clean up
	    for (int i=prevC.size()-2;i>=0;i--)
	    if (prevC[i])
	    {
	        if (fabs(b(Ts[i],Ts[i+1],T,q))<fr_eps)
		{
		    delete [] prevC[i];
		    prevC[i]=NULL;
		}
	    }
	}
	// corrector for dimensions for current T
	double corrector()
	{
	    if (q==1.0) return 1.0;
	    if (corrector_type==0)
		return 1;
	    if (corrector_type==1)
		return G2q*pow(T, q-1);
	    return 1;
	}
	// calculates rho for current T
	double rho(int i,int testing=0)
	{
//		double ka=47.91/(1.0+exp(106.06/(Tsolv->b_U[i]+(106.06/log(46.91))-T0))); // thermal coefficient
//		if (testing) ka=47.91/(1.0+exp(106.06/(Tsolv->testing_f(i*dL,T)+(106.06/log(46.91))-T0)));
		const double Q10=2.0;
		double ka=pow(Q10,0.1*(T0+273.15)*(Tsolv->b_U[i]-T0)/(Tsolv->b_U[i]+273.15)); // thermal coefficient
		if (testing) ka=pow(Q10,0.1*(T0+273.15)*(Tsolv->testing_f(i*dL,T)-T0)/(Tsolv->testing_f(i*dL,T)+273.15)); 
		double kb=fanox;
		double w10=0.5*(Hs->pr_vgm_s1s[i]+Hs->pr_vgm_s0s[i]);//Hs->wetness(i,10.1974);
		if (Hs->pr_w[i]<w10) 
		{
			kb=0.2+0.8*(Hs->pr_w[i]-Hs->pr_vgm_s0s[i])/(w10-Hs->pr_vgm_s0s[i]);
			if (testing)
			    kb=0.2+0.8*(Hs->wetness(0,-Hs->testing_f(i*dL,T))-Hs->pr_vgm_s0s[i])/(w10-Hs->pr_vgm_s0s[i]);
		}
		else
		{
			kb=1.0-(1.0-fanox)*(Hs->pr_w[i]-w10)/(Hs->pr_vgm_s1s[i]-w10);
			if (testing)
				kb=1.0-(1.0-fanox)*(Hs->wetness(0,-Hs->testing_f(i*dL,T))-w10)/(Hs->pr_vgm_s1s[i]-w10);
		}
		double kc=pr_Kc;
		return k[compound_id]*ka*kb*kc;
	}
	// out=rho*A*in for current T
	void Amult(double in[4],double out[4],int i,int testing=0)
	{
	    out[0]=-k[0]*in[0];
	    out[1]=-k[1]*in[1];
	    out[2]=alpha*(k[0]*in[0]+k[1]*in[1]+k[3]*in[3])+(alpha-1.0)*k[2]*in[2];
	    out[3]=beta*(k[0]*in[0]+k[1]*in[1]+k[2]*in[2])+(beta-1.0)*k[3]*in[3];
	    double r=rho(i,testing);
	    for (int j=0;j<4;j++)
		out[j]*=(r/pr_corr)*(k[j]/k[compound_id]);
	}
	// if CO2==1 - for CO2, if 0 - for CH4 (in C/s)
	double respiration(int CO2,int fout=0)
	{
	    double r1=0,r2=0,r,rch4; // r1=CO2, r2=CH4
	    double in[4],out[4];
	    for (int i=0;i<N+1;i++)
	    {
		//double lox=oxidation*exp(-ox_tau*i*dL); // oxidation in point
		in[compound_id]=b_U[i];
		for (int j=0;j<4;j++)
		    if (j!=compound_id)
			in[j]=Css[j]->b_U[i];
		Amult(in,out,i);
		// CO2 respiration in non-saturated state / CH4 in saturated
		r=dL*pr_corr*(out[0]+out[1]+out[2]+out[3])/(30.0*86400.0);
		rch4=0;
		if (!isfinite(r)) r=0.0;
		if (Hs->pr_w[i]<Hs->pr_vgm_s1s[i])
		    r1+=r;
		else
		{
		    rch4=(1.0-ox)*exp(-ox_tau*i*dL)*r;
		    r-=rch4;
		    r2+=rch4;
		    r1+=r;
		}
		if ((fout)&&((i%zoutstep)==0))
		{
		    if (CO2==0) fprintf(ch4_out,"%g ",-rch4/dL);
		    if (CO2==1) fprintf(co2_out,"%g ",-r/dL);
		}
	    }
	    if (CO2==0) // CH4
		return -r2;
	    return -r1;
	}
	// the value of integral in fractional integral/derivative's L1 approximation
	double b(double t0,double t1,double t,double q)
	{
	    double v1=0,v2=0;
	    if (q==1)
	    {
		if (t1==t) return 1;
		return 0;
	    }
	    if (t!=t1) v1=pow(t-t1,1.0-q);
	    if (t!=t0) v2=pow(t-t0,1.0-q);
	    return v2-v1;
	}
	// sets current SOC input and Kc
	void precalc_values()
	{
		pr_SoC=linear_from_pair(SoCT,SoCF,nSoC);
		pr_Kc=linear_from_pair(KcT,KcF,nKc);
		pr_corr=corrector();
		if (T0==1e6) // calculate average yearly temperature by monthly averaging of input data
		{
		    T0=0;
		    double st=Tsolv->T;
		    for (int i=0;i<12;i++)
		    {
			Tsolv->T=i*30*86400;
			T0+=Tsolv->linear_from_pair(Tsolv->TaT,Tsolv->TaF,Tsolv->nTa);
		    }
		    Tsolv->T=st;
		    T0/=12;
		    if (debug_level==1) printf("T0 %g\n",T0);
		}
	}
	// testing
	double testing_f(double z,double t) // solution for testing mode
	{
		return t*t*z*z/(30.0*86400*30.0*86400)+1+compound_id;
	}
	double testing_dfdt(double z,double t)
	{
		return  (Gamma(3.0)/Gamma(3-q))*pow(t,2-q)*z*z/(30.0*86400*30.0*86400);
	}
	double testing_dfdz(double z,double t)
	{
		return  2*z*t*t/(30.0*86400*30.0*86400);
	}
	double testing_d2fdz2(double z,double t)
	{
		return  2*t*t/(30.0*86400*30.0*86400);
	}
	double testing_F(double z,double t)
	{
		double F=0;
		F=testing_dfdt(z,t)+(-pr_soil_D[compound_id][0]*testing_d2fdz2(z,t)
		    +Hs->pr_vgm_k[0]*Hs->Kr(Hs->wetness(0,-Hs->testing_f(z,t)),0)*(Hs->testing_dfdz(z,t)-1)*testing_dfdz(z,t)
		    +(Hs->pr_vgm_k[0]*Hs->Kr(Hs->wetness(0,-Hs->testing_f(z,t)),0)*Hs->testing_d2fdz2(z,t)+
			    Hs->testing_dkdh(z,t)*Hs->testing_dfdz(z,t)*(Hs->testing_dfdz(z,t)-1))*testing_f(z,t))/pr_corr;
		double in[4],out[4],sc=compound_id;
		for (int i=0;i<4;i++)
		{
			compound_id=i;
			in[i]=testing_f(z,t);
		}
		compound_id=sc;
		Amult(in,out,z/dL,1);
		F-=out[compound_id]/(30.0*86400.0); 
		return F;
	}
	// coefficients of three-diagonal linear equations system
	double A_(int i,double *pickard_prevU)
	{
		return -((pr_soil_D[compound_id][i]/pr_corr)/(dL*dL))-Hs->b_V[i]/(2.0*dL*pr_corr);
	}
	double B(int i,double *pickard_prevU)
	{
		return -((pr_soil_D[compound_id][i]/pr_corr)/(dL*dL))+Hs->b_V[i]/(2.0*dL*pr_corr);
	}
	double R(int i,double *pickard_prevU)
	{
		return (2.0*(pr_soil_D[compound_id][i]/pr_corr)/(dL*dL))+(Hs->b_Vd[i]/pr_corr)+(1.0/tau)*(b(Ts[Ts.size()-1],T,T,q)/G2q);
	}
	// right part 
	double Rp(int ii,double *b_Uold,double *pickard_prev)
	{
		double ret=(1.0/tau)*(b(Ts[Ts.size()-1],T,T,q)/G2q)*b_Uold[ii];
		// SoC inputs
		double input_mult[4]={_gamma,1.0-_gamma,0,0};
		if (testing_mode!=1)
			ret+=pr_SoC*input_mult[compound_id]*inputs_divisor*0.5*(exp(-(ii-1)*dL*inputs_exponent)+exp(-ii*dL*inputs_exponent))/
				(((1/inputs_exponent)-exp(-inputs_exponent*L)/inputs_exponent)*pr_corr); // exponental decay of inputs with depth
		// SoC transpiration
		double in[4],out[4];
		in[compound_id]=pickard_prev[ii];
		for (int i=0;i<4;i++)
		    if (i!=compound_id)
			in[i]=Css[i]->b_U[ii];
		Amult(in,out,ii);
		ret+=out[compound_id]/(30*86400.0); // monthly rate to per second rate
		// fractional derivative's discretization
		double fr=0;
		if (q!=1.0)
		for (int i=prevC.size()-2;i>=0;i--)
		{
		    if ((prevC[i]==NULL)||(prevC[i+1]==NULL)) 
			break;
		    double v=b(Ts[i],Ts[i+1],T,q)*(prevC[i+1][ii]-prevC[i][ii])/(Ts[i+1]-Ts[i]);
		    if (fabs(b(Ts[i],Ts[i+1],T,q))<fr_eps) // finite memory
			break;
		    fr+=v;
		}
		fr/=G2q;
		ret-=fr;
		if (testing_mode==1)
			ret+=testing_F(ii*dL,T);
		return ret;
	}
	// linear solver
	// matrix multiplication of vector UU
	void mmult(double *UU)
	{
		mmult_main(UU);
#pragma omp single
		{
			// bottom boundary condition - flux=0
			if (testing_mode!=3) MM[N]=UU[N]+((pr_A[N]+pr_B[N])/(pr_R[N]+pr_B[N]*(2*dL*Hs->b_V[N]/pr_corr)/(pr_soil_D[compound_id][N]/pr_corr)))*UU[N-1];  // second order condition through "ghost point"
			// transparent flow 
			if (testing_mode==3) MM[N]=UU[N]+((pr_A[N]+pr_B[N])/(pr_R[N]))*UU[N-1];  // second order condition through "ghost point"
			// upper boundary condition - flux
			MM[0]=UU[0]+((pr_A[0]+pr_B[0])/(pr_R[0]-pr_A[0]*(2*dL*Hs->b_V[0]/pr_corr)/(pr_soil_D[compound_id][0]/pr_corr)))*UU[1];  // second order condition through "ghost point"
			if (testing_mode==1)
				MM[0]=UU[0];
		}
	}
	// precalculations for current time step
	void __precalc(double *pickard_prevU)
	{
		if (steady_state)
		{
		    printf("no steady state solver for C equation\n");
		    exit(0);
		}
#pragma omp single
		{ // precalculate non-linear coefficients using the solution of previous iteration
		    memcpy(sb_U,b_U,sizeof(double)*(N+2));
		    memcpy(b_U,pickard_prevU,sizeof(double)*(N+2));
		    if (q!=1.0)
		    if (Ts[Ts.size()-1]!=T-tau)
			save_to_storage();
		    precalc_values();
		    if (pr_A==NULL) pr_A=new double[N+2];
		    if (pr_B==NULL) pr_B=new double[N+2];
		    if (pr_R==NULL) pr_R=new double[N+2];
		    memcpy(b_U,sb_U,sizeof(double)*(N+2));
	    }
#pragma omp barrier
	    // precalc linear system values
		__precalc_main(pickard_prevU);
#pragma omp barrier
#pragma omp single
		{
			// bottom boundary condition - flux=0
			if (testing_mode!=3) RP[N]=Rp(N,b_U,pickard_prevU)/(pr_R[N]+pr_B[N]*(2*dL*Hs->b_V[N]/pr_corr)/(pr_soil_D[compound_id][N]/pr_corr));  // second order condition through "ghost point"
			// transparent flow 
			if (testing_mode==3) RP[N]=Rp(N,b_U,pickard_prevU)/pr_R[N];  // second order condition through "ghost point"
			// upper boundary condition - flux=part of SoC input
			{
				double input_mult[4]={_gamma,1.0-_gamma,0,0};
				RP[0]=(Rp(0,b_U,pickard_prevU)-pr_A[0]*(2*dL*(1.0/(pr_soil_D[compound_id][0]/pr_corr))*
					    pr_SoC*input_mult[compound_id]*(1.0-inputs_divisor)))
					/(pr_R[0]-pr_A[0]*(2*dL*Hs->b_V[0]/pr_corr)/(pr_soil_D[compound_id][0]/pr_corr));
			}
			if (testing_mode==1)
			{
				RP[N]=(Rp(N,b_U,pickard_prevU)-pr_B[N]*(2*dL*(1.0/(pr_soil_D[compound_id][0]/pr_corr))*((pr_soil_D[compound_id][0]/pr_corr)*testing_dfdz(N*dL,T)-
					Hs->pr_vgm_k[0]*Hs->Kr(Hs->wetness(0,-Hs->testing_f(N*dL,T)),0)*(Hs->testing_dfdz(N*dL,T)-1.0)*testing_f(N*dL,T)/pr_corr)))/
						(pr_R[N]+pr_B[N]*(2*dL*Hs->b_V[N]/pr_corr)/(pr_soil_D[compound_id][N]/pr_corr));  // second order condition through "ghost point"
				RP[0]=testing_f(0,T);
			}
		}
	}
	// constructor
	C_solver(double _L,double _C0,double _q,int _ct,int _ci,double _fa,double _ie,double _id, double _g,double _ox,double _ox_tau,H_solver *_Hs,T_solver *_Ts) : 
		ox(_ox),ox_tau(_ox_tau),C0(_C0),q(_q),corrector_type(_ct),compound_id(_ci),fanox(_fa),inputs_exponent(_ie),inputs_divisor(_id),_gamma(_g),Hs(_Hs),Tsolv(_Ts),basic_solver(_L)
	{
	    G2q=Gamma(2.0-q);
	    G1q=Gamma(1.0+q);
	    // initial conditions
	    for (int i = 0;i < N + 1;i++)
	    {
		b_U[i] = C0;
		if (testing_mode==1)
			b_U[i]=testing_f(i*dL,0);
	    }
	    save_to_storage();
	    char str[1024];
	    sprintf(str,"out_C%d.txt",compound_id);
	    out=fopen(str,"wt");
	    if (compound_id==3)
	    {
		rho_out=fopen("out_rho.txt","wt");
		soc_out=fopen("out_soc.txt","wt");
		co2_out=fopen("out_co2.txt","wt");
		ch4_out=fopen("out_ch4.txt","wt");
	    }
	    T0=1e6; // not set
	    Is.push_back(0);
	    R0s.push_back(0);
	    R1s.push_back(0);
	    Os.push_back(0);
	    fTs.push_back(0);
	    pr_corr=1;
	}
	// desctructor
	~C_solver()
	{
		if (SoCT) delete [] SoCT;
		if (SoCF) delete [] SoCF;
		for (int i=0;i<prevC.size();i++)
			delete [] prevC[i];
		fclose(out);
	        if (compound_id==3)
		{
			fclose(rho_out);
			fclose(soc_out);
			fclose(co2_out);
			fclose(ch4_out);
		}
	}
	// output
	void output_solution()
	{
		fprintf(out,"t(days) %g tau(seconds) %g SoC input %g C[%d]: ",(T-tau)/86400.0,tau,pr_SoC,compound_id);
		if (T==tau)
		{
		    for (int i=0;i<N+1;i+=zoutstep)
			fprintf(out,"%g ",i*dL);
		    fprintf(out,"\n");
		    fprintf(out,"t(days) %g tau(seconds) %g SoC input %g C[%d]: ",(T-tau)/86400.0,tau,pr_SoC,compound_id);
		}
		for (int i=0;i<N+1;i+=zoutstep)
		    fprintf(out,"%g ",b_U[i]);
		if (compound_id==3)
		if (Hs->pr_w)
		{
		    fprintf(rho_out,"t(days) %g rho: ",(T-tau)/86400.0);
		    fprintf(soc_out,"t(days) %g soc: ",(T-tau)/86400.0);
		    fprintf(co2_out,"t(days) %g co2_respiration: ",(T-tau)/86400.0);
		    fprintf(ch4_out,"t(days) %g ch4_respiration: ",(T-tau)/86400.0);
		    if (T==tau)
		    {
			for (int i=0;i<N+1;i+=zoutstep)
			{
			    fprintf(rho_out,"%g ",i*dL);
			    fprintf(soc_out,"%g ",i*dL);
			    fprintf(co2_out,"%g ",i*dL);
			    fprintf(ch4_out,"%g ",i*dL);
			}
			fprintf(rho_out,"\n");
			fprintf(soc_out,"\n");
			fprintf(co2_out,"\n");
			fprintf(ch4_out,"\n");
			fprintf(rho_out,"t(days) %g rho: ",(T-tau)/86400.0);
			fprintf(soc_out,"t(days) %g soc: ",(T-tau)/86400.0);
			fprintf(co2_out,"t(days) %g co2_respiration: ",(T-tau)/86400.0);
			fprintf(ch4_out,"t(days) %g ch4_respiration: ",(T-tau)/86400.0);
		    }
		    for (int i=0;i<N+1;i+=zoutstep)
		        fprintf(rho_out,"%g ",rho(i)/k[compound_id]);
		    fprintf(rho_out,"\n");
		    double r1=respiration(1,1);
		    double r2=respiration(0,1);
		    double sum=0;
		    double socn=0;
		    for (int i=0;i<N+1;i++)
		    {
			double lsum=dL*b_U[i];
			for (int j=0;j<4;j++)
			    if (j!=compound_id)
				lsum+=dL*Css[j]->b_U[i];
			sum+=lsum;
			if (i==N) socn=lsum/dL;
			if ((i%zoutstep)==0) fprintf(soc_out,"%g ",lsum/dL);
		    }
		    fprintf(soc_out,"\n");
		    fprintf(co2_out,"\n");
		    fprintf(ch4_out,"\n");
		    // balance
		    double If=0,R0f=0,R1f=0,Of=0;
		    if (T!=fTs[fTs.size()-1])
		    {
			Is.push_back(pr_SoC/pr_corr);
			R0s.push_back(r2/pr_corr);
			R1s.push_back(r1/pr_corr);
			Os.push_back((testing_mode==3)?Hs->b_V[N]*socn/pr_corr:0.0);
			fTs.push_back(T);
		    }
		    if (prev_out==0)
			balance=sum0=sum;
		    else
		    {
			// integrate fluxes
			for (int i=fTs.size()-2;i>=0;i--)
			{
			    double v,bv=b(fTs[i],fTs[i+1],T,1.0-q)/G1q;
			    v=bv*(Is[i+1]+Is[i])/2; If+=v;
			    v=bv*(R0s[i+1]+R0s[i])/2; R0f+=v;
			    v=bv*(R1s[i+1]+R1s[i])/2; R1f+=v;
			    v=bv*(Os[i+1]+Os[i])/2; Of+=v;
			    if (fabs(b(fTs[i],fTs[i+1],T,1.0-q))<fr_eps) // finite memory
				break;
			}
			balance=sum0+If-R0f-R1f-Of;
		    }
		    fprintf(out,"CO2 respiration - %g CH4 respiration - %g Balance %g sumSOC %g sumSoCbalance %g sumI %g sumRC02 %g sumRCH4 %g sumOut %g",r1,r2,pr_SoC-r1-r2,sum,balance,If,R1f,R0f,Of);
	        }
		if (testing_mode==1)
		{
		    double diff=0,sq1=0,sq2=0;
		    for (int i=0;i<N+1;i++)
		    {
			diff+=(b_U[i]-testing_f(i*dL,T-tau))*(b_U[i]-testing_f(i*dL,T-tau));
			sq1+=b_U[i]*b_U[i];
			sq2+=testing_f(i*dL,T-tau)*testing_f(i*dL,T-tau);
		    }
		    diff=sqrt(diff)/N;
		    sq1=sqrt(sq1);
		    sq2=sqrt(sq2);
		    printf ("testing diff %g rel %g\n",diff,fabs(sq1-sq2)/sq2);
		}
		fprintf(out,"\n");
		prev_out=T-tau;
	}
};
////////////////////////////////////////
// reads data from <time,value> pair file into rlT,rlV (nrl - their size) from filename
////////////////////////////////////////
void read_pair_file(double *&rlT,double *&rlV,int &nrl,char *filename)
{
    FILE *fi=fopen(filename,"rt");
    if (fi)
    {
	char str[1024];
	nrl=0;
	while (fgets(str,1024,fi)) nrl++;
	rlT=new double[nrl];
	rlV=new double[nrl];
	nrl=0;
	fseek(fi,0,SEEK_SET);			
	while (fscanf(fi,"%lg %lg\n",rlT+nrl,rlV+nrl)==2)
	    if (rlT[nrl]>=0)
		nrl++;
	printf("%d pairs read from %s\n",nrl,filename);
    }
    else
	printf("error reading %s\n",filename);
}
////////////////////////////////////////
////////////////// main/////////////////
////////////////////////////////////////
int main(int argc,char **argv)
{
	double Tm = 365*24*3600.0; // ending time, s
	double Om = 30*24*3600.0; // interval for solution output
	double H0=-1; // initial H, m water
	double L=3; // domain depth, m
	double _q=1.0;
	int _ct=0;
	double _fa=0.025;
	double _ie=5;
	double _id=0.7;
	double _g=0.4;
	double _ox=0.25;
	double _ox_tau=0.0;
	double C0[4]={1,2,3,4};
	// reading basic parameters
	for (int i=1;i<argc;i+=2)
	{
		// basic parameters
		if (strcmp(argv[i],"Tm")==0) // in days
			Tm=atof(argv[i+1])*86400;
		if (strcmp(argv[i],"Om")==0) // in days
			Om=atof(argv[i+1])*86400;
		if (strcmp(argv[i],"NB")==0)
			N=atoi(argv[i+1]);
		if (strcmp(argv[i],"H0")==0)
			H0=atof(argv[i+1]);
		if (strcmp(argv[i],"LB")==0)
			L=atof(argv[i+1]);
		if (strcmp(argv[i],"testing")==0)
			testing_mode=atoi(argv[i+1]);
		if (strcmp(argv[i],"debug_level")==0)
			debug_level=atoi(argv[i+1]);
		if (strcmp(argv[i],"Zoutstep")==0)
			zoutstep=atoi(argv[i+1]);
		// C equation parameters
		if (strcmp(argv[i],"C_q")==0)
			_q=atof(argv[i+1]);
		if (strcmp(argv[i],"corrector_type")==0)
			_ct=atoi(argv[i+1]);
		if (strcmp(argv[i],"fanox")==0)
			_fa=atof(argv[i+1]);
		if (strcmp(argv[i],"inputs_exponent")==0)
			_ie=atof(argv[i+1]);
		if (strcmp(argv[i],"inputs_divisor")==0)
			_id=atof(argv[i+1]);
		if (strcmp(argv[i],"oxidation")==0)
			_ox=atof(argv[i+1]);
		if (strcmp(argv[i],"ox_tau")==0)
			_ox_tau=atof(argv[i+1]);
		if (strcmp(argv[i],"gamma")==0)
			_g=atof(argv[i+1]);
		if (strcmp(argv[i],"C00")==0)
			C0[0]=atof(argv[i+1]);
		if (strcmp(argv[i],"C01")==0)
			C0[1]=atof(argv[i+1]);
		if (strcmp(argv[i],"C02")==0)
			C0[2]=atof(argv[i+1]);
		if (strcmp(argv[i],"C03")==0)
			C0[3]=atof(argv[i+1]);
		// fractional discretization parameters
		if (strstr(argv[i],"fr_eps")!=NULL)
		    fr_eps=atof(argv[i+1]);
		// linear solver parameters
		if (strstr(argv[i],"ls_eps")!=NULL)
		    ls_eps=atof(argv[i+1]);
		if (strstr(argv[i],"ls_max_iter")!=NULL)
		    ls_max_iter=atoi(argv[i+1]);
		if (strstr(argv[i],"ls_percent")!=NULL)
		    ls_percent=atof(argv[i+1]);
		if (strstr(argv[i],"ls_min_tau")!=NULL) //  in seconds
		    ls_min_tau=atof(argv[i+1]);
		if (strstr(argv[i],"ls_mult")!=NULL)
		    ls_mult=atof(argv[i+1]);
		if (strstr(argv[i],"max_tau")!=NULL) // in seconds
		    max_tau=atof(argv[i+1]);
	}
	printf("Grid size - %d, Tend %g Tsave %g \n",N,Tm,Om);
	// creating solvers
	H_solver *solver=new H_solver(L,H0);
	T_solver *Tsolver=new T_solver(L,solver);
	C_solver *Csolvers[4];
	for (int i=0;i<4;i++) Csolvers[i]=new C_solver(L,C0[i],_q,_ct,i,_fa,_ie,_id,_g,_ox,_ox_tau,solver,Tsolver);
	for (int i=0;i<4;i++) for (int j=0;j<4;j++) Csolvers[i]->Css[j]=Csolvers[j]; 
	// reading inputs and passing them to the solver
	int check=0;
	for (int i=1;i<argc;i+=2)
	{
		// reading the file with van Genuchtem-Mualem model's parameters
		if (strstr(argv[i],"vgm")!=NULL)
		{
		    FILE *fi=fopen(argv[i+1],"rt");
		    if (fi)
		    {
			char str[1024];
			vgm_nlayers=-1;
			while (fgets(str,1024,fi)) vgm_nlayers++;
			vgm_ns=new double[vgm_nlayers];
			vgm_as=new double[vgm_nlayers];
			vgm_s0s=new double[vgm_nlayers];
			vgm_s1s=new double[vgm_nlayers];
			vgm_h0=new double[vgm_nlayers];
			vgm_h1=new double[vgm_nlayers];
			vgm_k=new double[vgm_nlayers];
			vgm_power=new double[vgm_nlayers];
			vgm_specific_storage=new double[vgm_nlayers];
			soil_lambda=new double[vgm_nlayers];
			soil_Ct=new double[vgm_nlayers];
			for (int i=0;i<4;i++) soil_D[i]=new double[vgm_nlayers];
			vgm_nlayers=0;
			fseek(fi,0,SEEK_SET);
			fgets(str,1024,fi);
			while(fgets(str,1024,fi))
			{
				sscanf(str,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",vgm_s0s+vgm_nlayers,vgm_s1s+vgm_nlayers,
				vgm_ns+vgm_nlayers,vgm_as+vgm_nlayers,vgm_h0+vgm_nlayers,vgm_h1+vgm_nlayers,vgm_k+vgm_nlayers,vgm_power+vgm_nlayers,
				vgm_specific_storage+vgm_nlayers,soil_lambda+vgm_nlayers,soil_Ct+vgm_nlayers,
				soil_D[0]+vgm_nlayers,soil_D[1]+vgm_nlayers,soil_D[2]+vgm_nlayers,soil_D[3]+vgm_nlayers);
				vgm_nlayers++;
			}
			printf("VGM data read (%d layers)\n",vgm_nlayers);
			check++;
		    }
		}
		// ET
		if (strstr(argv[i],"et_file")!=NULL)
		{
			read_pair_file(solver->EV_T,solver->EV_F,solver->nEV_T,argv[i+1]);
			check++;
		}
		// precipitation
		if (strstr(argv[i],"prec_file")!=NULL)
		{
			read_pair_file(solver->perc_T,solver->perc_A,solver->nperc,argv[i+1]);
			check++;
		}
		// groundwater level
		if (strstr(argv[i],"gw_file")!=NULL)
		{
			read_pair_file(solver->flT,solver->flV,solver->nfl,argv[i+1]);
			check++;
		}
		// atmospheric temterature
		if (strstr(argv[i],"Ta_file")!=NULL)
		{
			read_pair_file(Tsolver->TaT,Tsolver->TaF,Tsolver->nTa,argv[i+1]);
			check++;
		}
		// SOC inputs
		if (strstr(argv[i],"SoC_file")!=NULL)
		{
			for (int ii=0;ii<4;ii++) read_pair_file(Csolvers[ii]->SoCT,Csolvers[ii]->SoCF,Csolvers[ii]->nSoC,argv[i+1]);
			check++;
		}
		// Kc
		if (strstr(argv[i],"Kc_file")!=NULL)
		{
			for (int ii=0;ii<4;ii++) read_pair_file(Csolvers[ii]->KcT,Csolvers[ii]->KcF,Csolvers[ii]->nKc,argv[i+1]);
			check++;
		}
	}
	if (check!=7)
	{
		printf("not all input file were given\n");
		return 1;
	}
	solver->soil_calc_coefs();
	Tsolver->soil_calc_coefs();
	Tsolver->initialize();
	for (int i=0;i<4;i++) Csolvers[i]->soil_calc_coefs();
	// solve steady state problem for H and T equations
	double smi=ls_max_iter;
	if (testing_mode!=1)
	{
	    solver->steady_state=1;
	    Tsolver->steady_state=1;	
	    ls_max_iter*=20;

	    solver->pickard_calc_step();
	    solver->calculate_velocity();
	    Tsolver->solve_up_to_T(solver->T-solver->tau,solver->tau);

	    solver->steady_state=0;
	    Tsolver->steady_state=0;	
	    solver->T=solver->tau;
	    Tsolver->T=Tsolver->tau;
	    ls_max_iter=smi;
	}
	// run simulation
	double last_o=0;
	solver->output_solution();
	Tsolver->output_solution();
	for (int i=0;i<4;i++) Csolvers[i]->output_solution();
	while (solver->T<Tm)
	{
		solver->pickard_calc_step();
		solver->calculate_velocity();
		Tsolver->solve_up_to_T(solver->T-solver->tau,solver->tau);
		Csolvers[3]->solve_up_to_T(solver->T-solver->tau,solver->tau,(basic_solver **)Csolvers,3);
		if ((solver->T-solver->tau-last_o)>=Om)
		{
		    last_o+=Om;
		    solver->output_solution();
		    Tsolver->output_solution();
		    for (int i=0;i<4;i++) Csolvers[i]->output_solution();
		    printf("T %g\n",solver->T-solver->tau);
		}
	}
	solver->output_solution();
	Tsolver->output_solution();
	for (int i=0;i<4;i++) Csolvers[i]->output_solution();
	return 0;
}
