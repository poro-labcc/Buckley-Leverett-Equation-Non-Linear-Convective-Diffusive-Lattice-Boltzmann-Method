//#include "pch.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <time.h>
#include <cmath>
#include <fstream>
#include <string>
using namespace std;
#pragma warning(disable:4996)
#define nc 10


int main()
{
	long int kk = 0, mn = 0, i = 0, j = 0, k = 0, l = 0, xx = 0, yy = 0, zz = 0;
	long int Nx = 201, mstep = 1000, Nx1=24300;
	int data[5], data1[nc][5], ncc=0, inc = 0, c;
	char mcname[1024];
	char velname[1024];
	char denname[1024];
	double w[5];
	int	cx[] = {0, 1, -1, 2, -2};
	double Asum,Umed,MassC,ncount;
	double time1,time2,mlups,Conv,feq,ncxsum,ncysum,nczsum,ncm;
	double Bx,t2A,GA_xx,PsiA,LaA,PiA_x,PiA_xx,PiA_xxx,Aeq,Aneq;
	double GA_xxx,La2A,skA;
	FILE* fp;
	FILE* MC;
	//=============================================double array==========================================================
	double** A = (double**)malloc(5 * sizeof(double));
	double** Ar = (double**)malloc(5 * sizeof(double));
	double* sk = (double*)malloc(Nx1 * sizeof(double)); 
	for (k = 0; k < 5; k++) {
		A[k] = (double*)malloc(Nx1 * sizeof(double));
		Ar[k] = (double*)malloc(Nx1 * sizeof(double));  
	}
	//----------------------------------------------------------------------------------------------------------------------
	ifstream theFile("Casos.txt");
	ncc = 0;
	while (theFile >> data[0] >> data[1] >> data[2]) {
		data1[ncc][0] = data[0];
		data1[ncc][1] = data[1];
		data1[ncc][2] = data[2];
		printf("Re=%d\t m=%d\t c=%d\n", data1[ncc][0], data1[ncc][1], data1[ncc][2]);
		ncc = ncc + 1;
	}
	// Pause and wait for user input
	std::cout << "Press Enter to continue...";
	std::cin.get(); // Wait for Enter key press
	for (inc = 0; inc < nc; inc++) {
		Nx = data1[inc][1];
		double Re = double(data1[inc][0]);
		c = data1[inc][2];
		double rel = pow(3, c); ;
		//----------------------------------------------------------------------------------------------------------------------
		sprintf(mcname, "MassCons-%ld-c%d", (Nx), data1[inc][2]);
		MC = fopen(mcname, "w");
		//Dados Adimensionais
		double uo = 1.0;
		double ho = 3.0;
		double nuo = uo*ho/Re;
		double epo=nuo*nuo*5.0;
		//Dados LBM
		double dx = ho / (Nx);
		double dt = dx / rel;
		double ue=uo/rel;
		double nu=nuo*dt/(dx*dx);
		double ep=epo/(dx*dx);
		double cs = 1.0 / sqrt(1.0);
		double tau = (nu / (cs * cs)) + (1. / 2.);
		double beta=-(ep/(tau*tau-tau+1.0/6.0));
		mstep=int(0.5/dt);
		// mstep=100;
		printf("Nx=%li\n" , Nx);
		printf("mstep=%li\n" , mstep);
		printf("ue=%e\t rel=%e\n", ue, rel);
		printf("dx=%e\t dt=%e\n", dx, dt);
		printf("nu=%e\t tau=%e\n", nu, tau);
		printf("epe=%e\t beta=%e\n", ep, beta);

		//-------------------------------------D1Q5-------------------------------------------------------
		w[0]=6.0/12.0;w[1]=2.0/12.0;w[2]=2.0/12.0;w[3]=1.0/12.0;w[4]=1.0/12.0;
		//-------------------------------------Init-------------------------------------------------------
		for (i = 0; i < Nx; i++) {
			sk[i] = 0.0;
			if ((i+0.5)*3.0/Nx>=0.75 and (i+0.5)*3.0/Nx<=2.25){
				sk[i] = 0.66;
			}
			Bx=ue*sk[i]*sk[i]/(sk[i]*sk[i]+(1.0-sk[i])*(1.0-sk[i]));
			for (k = 0; k < 5; k++) {
				A[k][i] = w[k]*(sk[i] + Bx*cx[k]);
			}
		}
		time1 = clock();
		//============================================MainLoop=============================================================
		for (kk = 1; kk <= mstep; kk++) {
			//--------------------------------------Colisao---------------------------------------------------
			for (i = 0; i < Nx; i++) {
				skA=sk[i];
				t2A=skA*skA+0.5*(1.0-skA)*(1.0-skA);
				Bx = ue*skA*skA/(t2A);
				GA_xx=ue*ue*( (486.0*skA*skA*skA*skA*skA-810.0*skA*skA*skA*skA+700.0*skA*skA*skA-276.0*skA*skA +102.0*skA -26.0)/(72.0*
						(27.0*skA*skA*skA*skA*skA*skA - 54.0*skA*skA*skA*skA*skA + 63.0*skA*skA*skA*skA - 44.0*skA*skA*skA +
						 21.0*skA*skA - 6.0*skA + 1.0)) + 27.0*sqrt(2.0)/72.0*atan(sqrt(2)*(3.0*skA-1.0)/2.0) );
				LaA=GA_xx+cs*cs*skA;
				GA_xxx= (ue*ue*ue*( (21870.0*pow(skA,9.0)-65610.0*pow(skA,8.0)+110160.0*pow(skA,7.0)-120960.0*pow(skA,6.0) + 97180.0*pow(skA,5.0)
                    -58660.0*pow(skA,4.0)+26560.0*pow(skA,3.0)-8240.0*pow(skA,2.0)+1670.0*skA-194.0)/(240.0*(243.0*pow(skA,10.0) 
                    -810.0*pow(skA,9.0)+1485.0*pow(skA,8.0)-1800.0*pow(skA,7.0)+1590.0*pow(skA,6.0)-1052.0*pow(skA,5.0)+503.0*pow(skA,4.0)
                    -200.0*pow(skA,3.0) + 55.0*pow(skA,2.0) - 10.0*skA + 1.0))
                    + 135.0*sqrt(2.0)/240.0*atan(sqrt(2.0)*(3.0*skA-1.0)/2.0) ) 
         			+ ue*2.0*(2.0*skA-1.0)/(3.0*skA*skA-2.0*skA + 1.0) );
				La2A=GA_xxx+cs*cs*cs*beta*(-Bx);
				PiA_x=0.0;
				PiA_xx=0.0;
				PiA_xxx=0.0;
				for (k = 0; k < 5; k++) {
					Aeq = w[k]*(skA + Bx*cx[k] +0.5*(LaA-skA)*(cx[k]*cx[k]-1.0) +0.5*(La2A-3.0*Bx)*(cx[k]*cx[k]-3.0)*cx[k] );
					PiA_x = PiA_x + (A[k][i]-Aeq)*cx[k];
					PiA_xx = PiA_xx + (A[k][i]-Aeq)*(cx[k]*cx[k]-1.0);
					PiA_xxx = PiA_xxx + (A[k][i]-Aeq)*cx[k]*(cx[k]*cx[k]-3.0);
				}
				for (k = 0; k < 5; k++) {
					Aeq = w[k]*(skA + Bx*cx[k] +0.5*(LaA-skA)*(cx[k]*cx[k]-1.0) +0.5*(La2A-3.0*Bx)*(cx[k]*cx[k]-3.0)*cx[k] );
					Ar[k][i] = Aeq + (1.0-1.0/tau)*w[k]*(cx[k]*PiA_x/(cs*cs));
				}
			}
			//---------------------------------------Streaming-FHH-BC-------------------------------------------------
			for (i = 0; i < Nx; i++) {
				for (k = 0; k < 5; k++) {
					xx = i + int(cx[k]);
					if (xx == Nx) { xx = 0; }
					if (xx == Nx+1) { xx = 1; }
					if (xx == -1) { xx = Nx-1; }
					if (xx == -2) { xx = Nx-2; }
					A[k][xx] = Ar[k][i];
				}
			}
			//**********************************Macro*************************************************
			for (i = 0; i < Nx; i++) {
				Asum = 0.0;
				for (k = 0; k < 5; k++) {
					Asum = Asum + A[k][i];
				}
				// sk[i]=Asum;
				sk[i]=Asum;
			}
			//----------------------------------Convergence------------------------------------------
			if (kk%10000==0 && mstep>=1000 || kk%10==0 && mstep<1000){
				ncount=0.0;
				MassC=0.0;
				for (i = 0; i < Nx; i++) {
					ncount=ncount+1.0;
					MassC = MassC + sk[i];
				}
				MassC = MassC/ncount;
				Umed=Umed/ncount;
				//Da=nu*Umed/(dp2*(Ny-2)*(Ny-2))*12.0;
				//Conv=Da-Da0;
				//Da0=Da;
				fprintf(MC, "%e ", MassC);
				time2 = clock();
				mlups = float(Nx * kk) / (1000000. * (time2 - time1) / 1000);
				//printf("kk=%li\t Umed=%f\t Da=%f\t Conv=%f\t mlups=%f\n", kk, Umed / ue, Da, Conv, mlups);
				printf("kk=%li\t MassC=%f\t mlups=%f\n", kk, MassC, mlups);
				//exit(0);
			}
			if (kk%mstep==0){
				ncount=0.0;
				MassC=0.0;
				for (i = 0; i < Nx; i++) {
					ncount=ncount+1.0;
					MassC = MassC + sk[i];
				}
				Umed=Umed/ncount;
				//Da=nu*Umed/(dp2*(Ny-2)*(Ny-2))*12.0;
				//Conv=Da-Da0;
				//Da0=Da;
				fprintf(MC, "%e ", MassC / ncount);
				time2 = clock();
				mlups = float(Nx * kk) / (1000000. * (time2 - time1) / 1000);
				//printf("kk=%li\t Umed=%f\t Da=%f\t Conv=%f\t mlups=%f\n", kk, Umed / ue, Da, Conv, mlups);
				printf("kk=%li\t MassC=%f\t mlups=%f\n", kk, MassC, mlups);

				sprintf(denname, "Densidade-%ld-c%d", (Nx), data1[inc][2]);
				fp = fopen(denname, "w");
				for (i = 0; i < Nx; i++) { fprintf(fp, "%e ", sk[i]); }
				fclose(fp);

				fclose(MC);
				//exit(0);
			}
		
		}
		sprintf(denname, "Densidade-%ld-c%d", (Nx), data1[inc][2]);
				fp = fopen(denname, "w");
				for (i = 0; i < Nx; i++) { fprintf(fp, "%e ", sk[i]); }
		fclose(fp);
	}
}

