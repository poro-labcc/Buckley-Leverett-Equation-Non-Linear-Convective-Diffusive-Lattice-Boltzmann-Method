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
#define nc 25


int main()
{
	long int kk = 0, mn = 0, i = 0, j = 0, k = 0, l = 0, xx = 0, yy = 0, zz = 0;
	long int Nx = 201, mstep = 1000, Nx1=40001;
	int data[3], data1[nc][3], ncc=0, inc = 0, c;
	char velname[1024];
	char denname[1024];
	double w[3];
	int	cx[] = {0, 1, -1};
	double Re,nuo,ho,uo,dx,rel,dt,nu,rhoin,cs,tau,omega,ue,dp,dp2,time1,t1,t2,force,Asum,Bsum,Vxsum,Vysum,Vzsum,Umed,Da0,Da,MassC,ncount;
	double time2,mlups,Conv,feq,Rr,beta,ncxsum,ncysum,nczsum,t3,ncm;
	double Bx, Bx2,t2A,t2B,GA_xx,GB_xx,PsiA,PsiB,LaA,LaB,PiA_x,PiB_x,Aeq,Aneq,Beq,Bneq;
	FILE* fp;
	FILE* MC;
	//=============================================double array==========================================================
	double** A = (double**)malloc(3 * sizeof(double));
	double** Ar = (double**)malloc(3 * sizeof(double));
	double** B = (double**)malloc(3 * sizeof(double));
	double** Br = (double**)malloc(3 * sizeof(double)); 
	double* sk = (double*)malloc(Nx1 * sizeof(double)); 
	double* sk2 = (double*)malloc(Nx1 * sizeof(double)); 
	for (k = 0; k < 3; k++) {
		A[k] = (double*)malloc(Nx1 * sizeof(double));
		Ar[k] = (double*)malloc(Nx1 * sizeof(double)); 
		B[k] = (double*)malloc(Nx1 * sizeof(double));
		Br[k] = (double*)malloc(Nx1 * sizeof(double)); 
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
		Re = double(data1[inc][0]);
		c = data1[inc][2];
		rel = pow(2, c); ;
		//----------------------------------------------------------------------------------------------------------------------
		MC = fopen("MassCons", "w");
		//Dados Adimensionais
		uo = 1.0;
		ho = 20.0;
		nuo = uo*ho/Re;
		//Dados LBM
		dx = ho / (Nx-1);
		dt = dx / rel;
		ue=uo/rel;
		nu=nuo*dt/(dx*dx);
		cs = 1.0 / sqrt(3.);
		tau = (nu / (cs * cs)) + (1. / 2.);
		mstep=int(10/dt);
		printf("Nx=%li\n" , Nx);
		printf("mstep=%li\n" , mstep);
		printf("ue=%e\t rel=%e\n", ue, rel);
		printf("dx=%e\t dt=%e\n", dx, dt);

		//-------------------------------------D1Q3-------------------------------------------------------
		w[0]=4.0/6.0;w[1]=1.0/6.0;w[2]=1.0/6.0;
		//-------------------------------------Init-------------------------------------------------------
		for (i = 0; i < Nx; i++) {
			sk[i] = 0.0;
			sk2[i] = (1.0-sk[i]);
			Bx=ue*sk[i]*sk[i]/(sk[i]*sk[i]+(1.0-sk[i])*(1.0-sk[i]));
			Bx2 = ue*sk2[i]*sk2[i]/(sk2[i]*sk2[i]+(1.0-sk2[i])*(1.0-sk2[i]));
			for (k = 0; k < 3; k++) {
				// A[k][i] = w[k]*(sk[i] + Bx*3.0*cx[k]);
				// B[k][i] = w[k]*(sk2[i] + Bx2*3.0*cx[k]);
				A[k][i] = w[k]*(sk[i] );
				B[k][i] = w[k]*(sk2[i]);
			}
		}
		time1 = clock();
		//============================================MainLoop=============================================================
		for (kk = 1; kk <= mstep; kk++) {
			//--------------------------------------Colisao---------------------------------------------------
			for (i = 0; i < Nx; i++) {
				Bx = ue*sk[i]*sk[i]/(sk[i]*sk[i]+(1.0-sk[i])*(1.0-sk[i]));
				Bx2 = ue*sk2[i]*sk2[i]/(sk2[i]*sk2[i]+(1.0-sk2[i])*(1.0-sk2[i]));
				for (k = 0; k < 3; k++) {
					Ar[k][i] = w[k]*(sk[i] + Bx*3.0*cx[k]);
					Br[k][i] = w[k]*(sk2[i] + Bx2*3.0*cx[k]);
				}
			}
			//---------------------------------------Streaming-FHH-BC-------------------------------------------------
			for (i = 0; i < Nx; i++) {
				for (k = 0; k < 3; k++) {
					xx = i + int(cx[k]);
					if (xx > Nx-1) { xx = 0; }
					if (xx < 0) { xx = Nx-1; }
					A[k][xx] = Ar[k][i];
					B[k][xx] = Br[k][i];
				}
			}
			A[1][0]= 1.0 - A[0][0]-A[2][0];
			A[0][Nx-1]=A[0][Nx-2];
			A[1][Nx-1]=A[1][Nx-2];
			A[2][Nx-1]=A[2][Nx-2];
			B[1][Nx-1]= 1.0 - B[0][Nx-1]-B[2][Nx-1];
			B[0][0]=B[0][1];
			B[1][0]=B[1][1];
			B[2][0]=B[2][1];
			//**********************************Macro*************************************************
			for (i = 0; i < Nx; i++) {
				Asum = 0.0;Bsum = 0.0; 
				for (k = 0; k < 3; k++) {
					Asum = Asum + A[k][i];
					Bsum = Bsum + B[k][i];
				}
				// sk[i]=Asum;
				// sk2[i]=Bsum;
				sk[i]=Asum;
				sk2[i]=(1.0-Asum);
			}
			//----------------------------------Convergence------------------------------------------
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
				for (i = 0; i < Nx; i++) { fprintf(fp, "%e ", sk2[i]); }
				fclose(fp);

				fclose(MC);
				//exit(0);
			}
		
		}
		sprintf(denname, "Densidade-%ld-c%d", (Nx), data1[inc][2]);
				fp = fopen(denname, "w");
				for (i = 0; i < Nx; i++) { fprintf(fp, "%e ", sk[i]); }
				for (i = 0; i < Nx; i++) { fprintf(fp, "%e ", sk2[i]); }
		fclose(fp);
	}
}

