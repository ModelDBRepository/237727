/************************************************************/
/*                                                          */
/*                                                          */
/*                                                          */
/*      SenMotifMD7.c                                       */
/*                                                          */
/*                                                          */
/*                                                          */
/*      2018/1/18                                           */
/************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <signal.h>

#define N_assm  8        // N_assm+1 number of cell assemblies f0-7 
#define N_T     19       // N_T+1 number of cells within a cell assemlby  
#define DT      1e-4     // discrete time 
#define delay   500      // delay between Nsen and Nmot 
#define cmPY    500e-12  // capacitance: P 
#define cmSB    243e-12  // BS
#define cmLB    115e-12  // BL
#define cmG     10e-12   // Glia(Astrocyte)
#define gmPY    25.0e-9  // conductance 
#define gmSB    9.7e-9   // 
#define gmLB    8.2e-9   // 
#define gmG     20.0e-9  // 
#define gGap    20.0e-9  // gap-junction
#define gAMPA   0.5e-9   // conductance for AMPA-receptors 
#define gGABA   0.7e-9   // conductance for GABA-receptors 
#define gGABAb  1e-9     // for GABAb receptor 
#define K1      9e4      //  
#define K2      1.2      //  
#define K3      180      //  
#define K4      34       //  
#define Kd      100      //  
#define nBS     4        // binding sites 

#define UPYact    -0.01   // action potential 
#define UPYres    -0.065  // reset membrane pot  
#define USBact    -0.01   // 
#define USBres    -0.07   //    
#define ULBact    -0.01   // 
#define ULBres    -0.07   //  
#define UGres     -0.07   //   

#define u_AMPA     0.0    // reversal potential 
#define u_GABA    -0.08   //  
#define u_GABAb   -0.095  //  

#define steep_PY   280.0  // steepness of sigmoid function 
#define thres_PY   -0.035 // threshold 
#define steep_PY2  280.0  //  
#define thres_PY2  -0.026 //  
#define steep_SB   300.0  //  
#define thres_SB   -0.035 //     
#define steep_LB   300.0  //  
#define thres_LB   -0.031 //  
#define steep_LB2   300.0 // 
#define thres_LB2  -0.031 //  

#define alph_AMPA    1.1e6   // open AMPA-channels 
#define beta_AMPA    190.0   // 
#define alph_GABA    5.0e6   // open GABA-channels 
#define beta_GABA    180.0   // 
#define Glut_c       0.001   // glutamate conc.   
#define GABA_cS      0.001   // GABA conc. released from BS    
#define GABA_cL      0.001   // GABA conc. released from BL

double m_G=1e9;              // GABA transport coeficient 
#define uG_trn       -0.07   // revesal pot. of glia

#define theta_inp0   0       // input-freature (as theta) "f0" 
#define theta_inp1   1       // 
#define theta_inp2   2       //  
#define theta_inp3   3       //  
#define theta_inp4   4       //   
#define theta_inp5   5       //  
#define theta_inp6   6       //  
#define theta_inp7   7       // 

#define int_inp0_0     0e-10      // input to f0 of Nsen 
#define int_inp0_1     3.6e-10    //  
#define int_inp0_2     int_inp0_0 // int_inp0_0 
#define int_inp0_3     int_inp0_0 // int_inp0_0
#define int_inp0_4     int_inp0_0 // int_inp0_0 
#define int_inp0_5     int_inp0_0 // int_inp0_0 
#define int_inp0_6     int_inp0_0 // int_inp0_0 
#define int_inp0_7     int_inp0_0 // int_inp0_0 
#define inp_prob       1.0        // input probability [0,1]

#define int_inpV2    2.5e-10    // inout to Nmot  
double delta_T=6e2;             // extrasynaptic GABAa receptors 

int t;                          // time development  
#define OUT       onset_0-10000 //  
#define PAT       "ON"          // display current network activity pattern 
#define INTV      99            // display interval 
#define PERIOD    OUT+30000     // time evolution period 
#define onset_0   30000         // sitmulus onset for Nsen 
#define period_0  10000         //  
#define onset_V2  1000          // for Nmot
#define period_V2 9999999       //  

#define wLPPV1		2.0	    	// weight P-P in Nsen 
#define wLPPV2		2.0	        // Nmot

void dfsGABA_ext(void);            // extrasynaptic GABA concentration  
double I_GABA[N_assm+2][N_T+2];    //  
double GABA_ext[N_assm+2][N_T+2];  // extrasynaptic GABA concentration 
double GABA_extw[N_assm+2][N_T+2]; // for working 
double gamma_=5.0;                  // decay
double GABA_c0=10E-07;             // basal GABA-concentration
double GABAamb_max=4E-06;          // max
double GABAamb_min=0E-07;	       // min
double rEXT[N_assm+2][N_T+2];      // fraction of open channels for extrasynaptic GABA 
double difrEXT(int,int);           // calculation of difference  

double I_MGB1[N_assm+2];           // input to current Nsen 
double uPY1[N_assm+2][N_T+2];      // membrane pot. of P  
double uSB1[N_assm+2][N_T+2];      // B 
double uLB1[N_assm+2][N_T+2];      // BL 
double uG[N_assm+2][N_T+2];        // glia 
double vPY1[N_assm+2][N_T+2][PERIOD+2];  // action potential of P
double vSB1[N_assm+2][N_T+2][PERIOD+2];  // BS 
double vLB1[N_assm+2][N_T+2][PERIOD+2];  // BL 
double rPY1[N_assm+2][N_T+2];            // fraction of open channels for AMPA-receptors in P 
double rPY1d[N_assm+2][N_T+2][PERIOD+2]; // delay to Nmot 
double rSB1[N_assm+2][N_T+2];            // BS 
double sF[N_assm+2][N_T+2];              // rBS  
double rLB1[N_assm+2][N_T+2];            // GABA-receptors by BL
double GlutPY_c1[N_assm+2][N_T+2];       // glugamate release to P
double GlutSB_c1[N_assm+2][N_T+2];       // glugamate release to BS
double GlutLB_c1[N_assm+2][N_T+2];       // glugamate release to BL 
double duPY_leak1[N_assm+2][N_T+2];      // leak current
double duPY_rec_1[N_assm+2][N_T+2];      // recurrent excitatory current
double duPY_fed_1[N_assm+2][N_T+2];      // feedback inhibitory current
double duPY_lat_1[N_assm+2][N_T+2];      // lateral inhibitory current
double duPY_topdown[N_assm+2][N_T+2];    // Nmot-to-Nsen 
double duPY_ext_1[N_assm+2][N_T+2];      // tonic nhibitory current
double duPY_MGB1[N_assm+2];              // external current
double duSB_leak1[N_assm+2][N_T+2];     
double duSB_1[N_assm+2][N_T+2];
double duSB_ext_1;                       
double duLB_leak1[N_assm+2][N_T+2];      
double duLB1[N_assm+2][N_T+2];
double duG_leak[N_assm+2][N_T+2];        
double duG_Ia[N_assm+2][N_T+2];          // inhibitory current from BS
double duG_G[N_assm+2][N_T+2];           // gapjunction 
double drPY1[N_assm+2][N_T+2];           // fraction of open channels for AMPA-receptors by P 
double drSB1[N_assm+2][N_T+2];           // BS 
double drLB1[N_assm+2][N_T+2];           // BL 

double AMPA_c1[N_assm+2][N_T+2];                   // glutamate synpaptic concentration 
double GABASB_c1[N_assm+2][N_T+2][PERIOD+2];       // GABA synpaptic concentration from BS
double GABALB_c1[N_assm+2][N_T+2][PERIOD+2];       // BL  
double w_rec_1[N_assm+2][N_T+2][N_assm+2][N_T+2];  // P-to-P weight
double w_lat_1[N_assm+2][N_T+2][N_T+2];            // P-to-BL weight
double w_v1Ib_v2[N_assm+2][N_T+2][N_assm+2][N_T+2];// topdown V2(P)-to-V1(Ib) 結合荷重 
double w_v1P_v2[N_assm+2][N_T+2][N_assm+2][N_T+2]; // topdown V2(P)-to-V1(P) 結合荷重 
double w_v1G_v2[N_assm+2][N_T+2][N_assm+2][N_T+2]; // topdown V2(P)-to-V1(G) 結合荷重
double wSB_PY1[N_assm+2][N_T+2];                   // セルアセンブリ内局所抑制(SBC-to-PYC)結合荷重
double wLB_PY1[N_assm+2][N_assm+2][N_T+2];         // P-to-BL weight
double wG_Ia[N_assm+2];                            // BS-to-glia weight
double difuPY1(int,int);   // difference P membrane pot  
double difuSB1(int,int);   // BS  
double difuLB1(int,int);   // BL  
double difuG(int,int);     // glia  
double difrPY1(int,int);   // fraction of open channes for P
double difrSB1(int,int);   // BS  
double difsF(int,int);     // for rSB1
double difrLB1(int,int);   // BL  

//Nmot
double I_V2[N_assm+2];                                
double uPY2[N_assm+2][N_T+2];                        
double uLB2[N_assm+2][N_T+2];                      
double vPY2[N_assm+2][N_T+2][PERIOD+2];           
double vLB2[N_assm+2][N_T+2][PERIOD+2];        
double rPY2[N_assm+2][N_T+2];                    
double rPY2d[N_assm+2][N_T+2][PERIOD+2]; 
double rLB2[N_assm+2][N_T+2];           
double GlutPY_c2[N_assm+2][N_T+2];         
double GlutLB_c2[N_assm+2][N_T+2];                
double duPY_leak2[N_assm+2][N_T+2];             
double duPY_rec_2[N_assm+2][N_T+2];              
double duPY_fed_2[N_assm+2][N_T+2];              
double duPY_lat_2[N_assm+2][N_T+2];               
double duPY_bottomup[N_assm+2][N_T+2];            
double duPY_ext_2;                                 
double duPY_MGB2[N_assm+2];                        
double duPY_inp;                                 
double duLB_leak2[N_assm+2][N_T+2];                
double duLB2[N_assm+2][N_T+2];
double duLB2_bottom_up[N_assm+2][N_T+2];           
double duLB_ext2;                                
double drPY2[N_assm+2][N_T+2];                    
double drLB2[N_assm+2][N_T+2];                    
double AMPA_c2[N_assm+2][N_T+2];                   
double GABALB_c2[N_assm+2][N_T+2][PERIOD+2];      
double w_rec_2[N_assm+2][N_T+2][N_assm+2][N_T+2];  
double w_lat_2[N_assm+2][N_T+2][N_T+2];            
double w_v2_v1[N_assm+2][N_T+2][N_assm+2][N_T+2] ; 
double wSB_PY2[N_assm+2][N_T+2];                  
double wLB_PY2[N_assm+2][N_assm+2][N_T+2];         
double difuPY2(int,int);  
double difuLB2(int,int);    
double difrPY2(int,int);  
double difrLB2(int,int);   
double sigmoidPY(double);  // sigmoid for P of Nsen
double sigmoidPY2(double); // Nmot
double sigmoidSB(double); 
double sigmoidLB(double); 
double sigmoidLB2(double); 
double rand01(long int *); // radom number [0, 1.0]
void init(void);       
long int SEEDMP=10000; 

// Nsen
FILE *gaba0_01,*gaba1_01,*gaba2_01,*gaba3_01,*gaba4_01,*gaba5_01,*gaba6_01,*gaba7_01;
FILE *uPY0_01,*uPY1_01,*uPY2_01,*uPY3_01,*uPY4_01,*uPY5_01,*uPY6_01,*uPY7_01;
FILE *uSB0_01,*uSB1_01,*uSB2_01,*uSB3_01,*uSB4_01,*uSB5_01,*uSB6_01,*uSB7_01;
FILE *uLB0_01,*uLB1_01,*uLB2_01,*uLB3_01,*uLB4_01,*uLB5_01,*uLB6_01,*uLB7_01;
FILE *uG0_01,*uG1_01,*uG2_01,*uG3_01,*uG4_01,*uG5_01,*uG6_01,*uG7_01;

// Nmot
FILE *uPY0_02,*uPY1_02,*uPY2_02,*uPY3_02,*uPY4_02,*uPY5_02,*uPY6_02,*uPY7_02;
FILE *uLB0_02,*uLB1_02,*uLB2_02,*uLB3_02,*uLB4_02,*uLB5_02,*uLB6_02,*uLB7_02;


void main(void){
	int theta,i,ii;
	double sigPYT,sigSBT,sigLBT;
	int ActPY1_ON[N_assm+2][N_T+2],ActSB1_ON[N_assm+2][N_T+2],ActLB1_ON[N_assm+2][N_T+2],duration;
	int ActPY2_ON[N_assm+2][N_T+2],ActLB2_ON[N_assm+2][N_T+2];
	int ActPY1_OF[N_assm+2][N_T+2],ActSB1_OF[N_assm+2][N_T+2],ActLB1_OF[N_assm+2][N_T+2];
	int ActPY2_OF[N_assm+2][N_T+2],ActLB2_OF[N_assm+2][N_T+2];

	init(); 

	gaba0_01=fopen("gaba0_01.dat","w");
	gaba1_01=fopen("gaba1_01.dat","w");
	gaba2_01=fopen("gaba2_01.dat","w");
	gaba3_01=fopen("gaba3_01.dat","w");
	gaba4_01=fopen("gaba4_01.dat","w");
	gaba5_01=fopen("gaba5_01.dat","w");
	gaba6_01=fopen("gaba6_01.dat","w");
	gaba7_01=fopen("gaba7_01.dat","w");
    
	// Nsen
	uPY0_01=fopen("uPY0_01.dat","w");
	uPY1_01=fopen("uPY1_01.dat","w");
	uPY2_01=fopen("uPY2_01.dat","w");
	uPY3_01=fopen("uPY3_01.dat","w");
	uPY4_01=fopen("uPY4_01.dat","w");
	uPY5_01=fopen("uPY5_01.dat","w");
	uPY6_01=fopen("uPY6_01.dat","w");
	uPY7_01=fopen("uPY7_01.dat","w");

	uSB0_01=fopen("uSB0_01.dat","w");
	uSB1_01=fopen("uSB1_01.dat","w");
	uSB2_01=fopen("uSB2_01.dat","w");
	uSB3_01=fopen("uSB3_01.dat","w");
	uSB4_01=fopen("uSB4_01.dat","w");
	uSB5_01=fopen("uSB5_01.dat","w");
	uSB6_01=fopen("uSB6_01.dat","w");
	uSB7_01=fopen("uSB7_01.dat","w");

	uLB0_01=fopen("uLB0_01.dat","w");
	uLB1_01=fopen("uLB1_01.dat","w");
	uLB2_01=fopen("uLB2_01.dat","w");
	uLB3_01=fopen("uLB3_01.dat","w");
	uLB4_01=fopen("uLB4_01.dat","w");
	uLB5_01=fopen("uLB5_01.dat","w");
	uLB6_01=fopen("uLB6_01.dat","w");
	uLB7_01=fopen("uLB7_01.dat","w");

	uG0_01=fopen("uG0_01.dat","w");
	uG1_01=fopen("uG1_01.dat","w");
	uG2_01=fopen("uG2_01.dat","w");
	uG3_01=fopen("uG3_01.dat","w");
	uG4_01=fopen("uG4_01.dat","w");
	uG5_01=fopen("uG5_01.dat","w");
	uG6_01=fopen("uG6_01.dat","w");
	uG7_01=fopen("uG7_01.dat","w");

	// Nmot
	uPY0_02=fopen("uPY0_02.dat","w");
	uPY1_02=fopen("uPY1_02.dat","w");
	uPY2_02=fopen("uPY2_02.dat","w");
	uPY3_02=fopen("uPY3_02.dat","w");
	uPY4_02=fopen("uPY4_02.dat","w");
	uPY5_02=fopen("uPY5_02.dat","w");
	uPY6_02=fopen("uPY6_02.dat","w");
	uPY7_02=fopen("uPY7_02.dat","w");

	uLB0_02=fopen("uLB0_02.dat","w");
	uLB1_02=fopen("uLB1_02.dat","w");
	uLB2_02=fopen("uLB2_02.dat","w");
	uLB3_02=fopen("uLB3_02.dat","w");
	uLB4_02=fopen("uLB4_02.dat","w");
	uLB5_02=fopen("uLB5_02.dat","w");
	uLB6_02=fopen("uLB6_02.dat","w");
	uLB7_02=fopen("uLB7_02.dat","w");

	for (theta=0; theta<=N_assm; ++theta){  
		for (ii=0; ii<=N_T; ++ii){
			ActPY1_ON[theta][ii] = 0; 
			ActSB1_ON[theta][ii] = 0; 
			ActLB1_ON[theta][ii] = 0; 
			ActPY1_OF[theta][ii] = 0; 
			ActSB1_OF[theta][ii] = 0; 
			ActLB1_OF[theta][ii] = 0; 
			uPY1[theta][ii]=UPYres;
			uSB1[theta][ii]=USBres;
			uLB1[theta][ii]=ULBres;
		}
	}
	for (theta=0; theta<=N_assm; ++theta){  
		for (ii=0; ii<=N_T; ++ii){
			ActPY2_ON[theta][ii] = 0; 
			ActLB2_ON[theta][ii] = 0; 
			ActPY2_OF[theta][ii] = 0; 
			ActLB2_OF[theta][ii] = 0; 
			uPY2[theta][ii]=UPYres;
			uLB2[theta][ii]=ULBres;
		}
	}
	duration  = (int)(0.001/DT);
	for (t=0; t<PERIOD; ++t){    
		for (theta=0; theta<=N_assm-1; ++theta){
			for (i=0; i<=N_T; ++i){	
				uPY1[theta][i] += difuPY1(theta,i); 
				sigPYT = sigmoidPY(uPY1[theta][i]);
				if (ActPY1_ON[theta][i]==0 && ActPY1_OF[theta][i]==0) vPY1[theta][i][t]=(double)(int)(rand01(&SEEDMP)+sigPYT);
				if (vPY1[theta][i][t]==1.0){
					uPY1[theta][i]=UPYact;
					GlutPY_c1[theta][i]=Glut_c;
					ActPY1_ON[theta][i]=1;
					ActPY1_OF[theta][i]=1;
				}else{
					GlutPY_c1[theta][i]=0.0;
				}
				if (ActPY1_ON[theta][i] != 0){
					++ActPY1_ON[theta][i];
					uPY1[theta][i]=UPYact;
					GlutPY_c1[theta][i]=Glut_c;
				}
				if (ActPY1_OF[theta][i] != 0){
					++ActPY1_OF[theta][i];
				}
				if (ActPY1_ON[theta][i] > duration){
					uPY1[theta][i]=UPYres;
					ActPY1_ON[theta][i] = 0;
				}
				if (ActPY1_OF[theta][i] > 10*duration){
					ActPY1_OF[theta][i] = 0;
				}
				rPY1[theta][i] += difrPY1(theta,i);   
				rPY1d[theta][i][t] = rPY1[theta][i]; 
				if (t%INTV == 0){
				printf("time = %d\n",t);
				}
				uSB1[theta][i] += difuSB1(theta,i); 
				sigSBT = sigmoidSB(uSB1[theta][i]);
				if (ActSB1_ON[theta][i]==0 && ActSB1_OF[theta][i]==0) vSB1[theta][i][t]=(double)(int)(rand01(&SEEDMP)+sigSBT);
				if (vSB1[theta][i][t]==1.0){
					uSB1[theta][i]=USBact;
					GABASB_c1[theta][i][t]=GABA_cS;
					ActSB1_ON[theta][i]=1;
					ActSB1_OF[theta][i]=1;
				}else{
					GABASB_c1[theta][i][t]=0.0;
				}
				if (ActSB1_ON[theta][i] != 0){
					++ActSB1_ON[theta][i];
					uSB1[theta][i]=USBact;
					GABASB_c1[theta][i][t]=GABA_cS;
				}
				if (ActSB1_OF[theta][i] != 0){
					++ActSB1_OF[theta][i];
				}
				if (ActSB1_ON[theta][i] > duration){
					uSB1[theta][i]=USBres;
					ActSB1_ON[theta][i] = 0;
				}
				if (ActSB1_OF[theta][i] > 10*duration){
					ActSB1_OF[theta][i] = 0;
				}	
				rSB1[theta][i] += difrSB1(theta,i); 
				sF[theta][i]   += difsF(theta,i);    
				uLB1[theta][i] += difuLB1(theta,i); 
				sigLBT = sigmoidLB(uLB1[theta][i]);
				if (ActLB1_ON[theta][i]==0 && ActLB1_OF[theta][i]==0) vLB1[theta][i][t]=(double)(int)(rand01(&SEEDMP)+sigLBT);
				if (vLB1[theta][i][t]==1.0){
					uLB1[theta][i]=ULBact;
					GABALB_c1[theta][i][t]=GABA_cL;
					ActLB1_ON[theta][i]=1;
					ActLB1_OF[theta][i]=1;
				}else{
					GABALB_c1[theta][i][t]=0.0;
				}
				if (ActLB1_ON[theta][i] != 0){
					++ActLB1_ON[theta][i];
					uLB1[theta][i]=ULBact;
					GABALB_c1[theta][i][t]=GABA_cL;
				}
				if (ActLB1_OF[theta][i] != 0){
					++ActLB1_OF[theta][i];
				}				
				if (ActLB1_ON[theta][i] > duration){
					uLB1[theta][i]=ULBres;
					ActLB1_ON[theta][i] = 0;
				}
				if (ActLB1_OF[theta][i] > 10*duration){
					ActLB1_OF[theta][i] = 0;
				}
				rLB1[theta][i] += difrLB1(theta,i); 
				uG[theta][i] += difuG(theta,i); 
				uPY2[theta][i] += difuPY2(theta,i); 
				sigPYT = sigmoidPY2(uPY2[theta][i]);
				if (ActPY2_ON[theta][i]==0 && ActPY2_OF[theta][i]==0) vPY2[theta][i][t]=(double)(int)(rand01(&SEEDMP)+sigPYT);
				if (vPY2[theta][i][t]==1.0){
					uPY2[theta][i]=UPYact;
					GlutPY_c2[theta][i]=Glut_c;
					ActPY2_ON[theta][i]=1;
					ActPY2_OF[theta][i]=1;
				}else{
					GlutPY_c2[theta][i]=0.0;
				}
				if (ActPY2_ON[theta][i] != 0){
					++ActPY2_ON[theta][i];
					uPY2[theta][i]=UPYact;
					GlutPY_c2[theta][i]=Glut_c;
				}
				if (ActPY2_OF[theta][i] != 0){
					++ActPY2_OF[theta][i];
				}
				if (ActPY2_ON[theta][i] > duration){
					uPY2[theta][i]=UPYres;
					ActPY2_ON[theta][i] = 0;
				}
				if (ActPY2_OF[theta][i] > 10*duration){
					ActPY2_OF[theta][i] = 0;
				}
				rPY2[theta][i] += difrPY2(theta,i); 
				rPY2d[theta][i][t] = rPY2[theta][i]; 
				if (t%INTV == 0){
				}
				uLB2[theta][i] += difuLB2(theta,i); 
				sigLBT = sigmoidLB2(uLB2[theta][i]);
				if (ActLB2_ON[theta][i]==0 && ActLB2_OF[theta][i]==0) vLB2[theta][i][t]=(double)(int)(rand01(&SEEDMP)+sigLBT);
				if (vLB2[theta][i][t]==1.0){
					uLB2[theta][i]=ULBact;
					GABALB_c2[theta][i][t]=GABA_cL;
					ActLB2_ON[theta][i]=1;
					ActLB2_OF[theta][i]=1;
				}else{
					GABALB_c2[theta][i][t]=0.0;
				}
				if (ActLB2_ON[theta][i] != 0){
					++ActLB2_ON[theta][i];
					uLB2[theta][i]=ULBact;
					GABALB_c2[theta][i][t]=GABA_cL;
				}
				if (ActLB2_OF[theta][i] != 0){
					++ActLB2_OF[theta][i];
				}				
				if (ActLB2_ON[theta][i] > duration){
					uLB2[theta][i]=ULBres;
					ActLB2_ON[theta][i] = 0;
				}
				if (ActLB2_OF[theta][i] > 10*duration){
					ActLB2_OF[theta][i] = 0;
				}
				rLB2[theta][i] += difrLB2(theta,i); 

				if (t>=OUT){  
					if (theta==0 && i==0){									
						fprintf(uPY0_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY0_02,"%f\n",uPY2[theta][i]); 
						fprintf(uSB0_01,"%f\n",uSB1[theta][i]); 
						fprintf(uLB0_01,"%f\n",uLB1[theta][i]); 
						fprintf(uG0_01,"%f\n",uG[theta][i]); 
						fprintf(gaba0_01,"%e\n",GABA_ext[theta][i]);
					} 
					if (theta==1 && i==0){									
						fprintf(uPY1_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY1_02,"%f\n",uPY2[theta][i]); 
						fprintf(uSB1_01,"%f\n",uSB1[theta][i]); 
						fprintf(uLB1_01,"%f\n",uLB1[theta][i]); 
						fprintf(uLB1_02,"%f\n",uLB2[theta][i]); 
						fprintf(uG1_01,"%f\n",uG[theta][i]); 
						fprintf(gaba1_01,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==2 && i==0){									
						fprintf(uPY2_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY2_02,"%f\n",uPY2[theta][i]); 
						fprintf(uSB2_01,"%f\n",uSB1[theta][i]); 
						fprintf(uLB2_01,"%f\n",uLB1[theta][i]); 
						fprintf(uLB2_02,"%f\n",uLB2[theta][i]); 
						fprintf(uG2_01,"%f\n",uG[theta][i]); 
						fprintf(gaba2_01,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==3 && i==0){									
						fprintf(uPY3_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY3_02,"%f\n",uPY2[theta][i]); 
						fprintf(uSB3_01,"%f\n",uSB1[theta][i]); 
						fprintf(uLB3_01,"%f\n",uLB1[theta][i]); 
						fprintf(uLB3_02,"%f\n",uLB2[theta][i]); 
						fprintf(uG3_01,"%f\n",uG[theta][i]); 
						fprintf(gaba3_01,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==4 && i==0){									
						fprintf(uPY4_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY4_02,"%f\n",uPY2[theta][i]); 
						fprintf(uSB4_01,"%f\n",uSB1[theta][i]); 
						fprintf(uLB4_01,"%f\n",uLB1[theta][i]); 
						fprintf(uLB4_02,"%f\n",uLB2[theta][i]); 
						fprintf(uG4_01,"%f\n",uG[theta][i]); 
						fprintf(gaba4_01,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==5 && i==0){									
						fprintf(uPY5_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY5_02,"%f\n",uPY2[theta][i]); 
						fprintf(uSB5_01,"%f\n",uSB1[theta][i]); 
						fprintf(uLB5_01,"%f\n",uLB1[theta][i]); 
						fprintf(uG5_01,"%f\n",uG[theta][i]); 
						fprintf(gaba5_01,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==6 && i==0){									
						fprintf(uPY6_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY6_02,"%f\n",uPY2[theta][i]); 
						fprintf(uSB6_01,"%f\n",uSB1[theta][i]); 
						fprintf(uLB6_01,"%f\n",uLB1[theta][i]); 
						fprintf(uG6_01,"%f\n",uG[theta][i]); 
						fprintf(gaba6_01,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==7 && i==0){									
						fprintf(uPY7_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY7_02,"%f\n",uPY2[theta][i]); 
						fprintf(uSB7_01,"%f\n",uSB1[theta][i]); 
						fprintf(uLB7_01,"%f\n",uLB1[theta][i]); 
						fprintf(uG7_01,"%f\n",uG[theta][i]); 
						fprintf(gaba7_01,"%e\n",GABA_ext[theta][i]);
					}
				} 
			} 
		} 

		dfsGABA_ext(); 
		for (theta=0; theta<=N_assm; ++theta){
			for (i=0; i<=N_T; ++i){	
			rEXT[theta][i] += difrEXT(theta,i); 
			}
		}
	} 
	fclose(gaba0_01);
	fclose(gaba1_01);
	fclose(gaba2_01);
	fclose(gaba3_01);
	fclose(gaba4_01);
	fclose(gaba5_01);
	fclose(gaba6_01);
	fclose(gaba7_01);
    
	fclose(uPY0_01);
	fclose(uPY1_01);
	fclose(uPY2_01);
	fclose(uPY3_01);
	fclose(uPY4_01);
	fclose(uPY5_01);
	fclose(uPY6_01);
	fclose(uPY7_01);

	fclose(uSB0_01);
	fclose(uSB1_01);
	fclose(uSB2_01);
	fclose(uSB3_01);
	fclose(uSB4_01);
	fclose(uSB5_01);
	fclose(uSB6_01);
	fclose(uSB7_01);

	fclose(uLB0_01);
	fclose(uLB1_01);
	fclose(uLB2_01);
	fclose(uLB3_01);
	fclose(uLB4_01);
	fclose(uLB5_01);
	fclose(uLB6_01);
	fclose(uLB7_01);

	fclose(uG0_01);
	fclose(uG1_01);
	fclose(uG2_01);
	fclose(uG3_01);
	fclose(uG4_01);
	fclose(uG5_01);
	fclose(uG6_01);
	fclose(uG7_01);

	fclose(uPY0_02);
	fclose(uPY1_02);
	fclose(uPY2_02);
	fclose(uPY3_02);
	fclose(uPY4_02);
	fclose(uPY5_02);
	fclose(uPY6_02);
	fclose(uPY7_02);

	fclose(uLB0_02);
	fclose(uLB1_02);
	fclose(uLB2_02);
	fclose(uLB3_02);
	fclose(uLB4_02);
	fclose(uLB5_02);
	fclose(uLB6_02);
	fclose(uLB7_02);

	printf("\a");
    printf("\a");
    printf("\a");
} 

double difuPY1(int thetaa, int ii){
	int jj,thetdash;
	double duPY1=0.0,rec_exT1=0.0,lat_exT=0.0,lat_ihTT=0.0,lat_ihTP=0.0,exP=0.0,c0=0.5,c1=1.0,pw=2.0;
	double top_down = 0.0;
	double alpha_T=1.0;  
	double tau_T0=0.1,tau_T1=0.1,tau_T2=0.01,tau_T3=0.01;    

	duPY_leak1[thetaa][ii] = -(DT*gmPY/cmPY)*(uPY1[thetaa][ii]-UPYres); 
	for (jj=0; jj<N_T; ++jj){ 
		rec_exT1 += w_rec_1[thetaa][ii][thetaa][jj]*rPY1[thetaa][jj]; 
	}

	for (thetdash=0; thetdash<=N_assm-1; ++thetdash){
		for (jj=0; jj<=N_T; ++jj){
			if (thetaa!=thetdash) rec_exT1 += w_rec_1[thetaa][ii][thetdash][jj]*rPY1[thetdash][jj];
		}
	}

	duPY_rec_1[thetaa][ii] = -(DT/cmPY)*gAMPA*(uPY1[thetaa][ii]-u_AMPA)*rec_exT1; 
	for (jj=0; jj<N_T; ++jj){ 
		lat_ihTT += w_lat_1[thetaa][ii][jj]*rLB1[thetaa][jj];
	}

	duPY_lat_1[thetaa][ii] = -(DT/cmPY)*gGABA*(uPY1[thetaa][ii]-u_GABA)*lat_ihTT; 
	for (thetdash=0; thetdash<=N_assm-1; ++thetdash){
		for (jj=0; jj<=N_T; ++jj){
			top_down += w_v1P_v2[thetaa][ii][thetdash][jj]*rPY2d[thetdash][jj][t-delay];
		}
	}

	duPY_topdown[thetaa][ii] = -(DT/cmPY)*gAMPA*(uPY1[thetaa][ii]-u_AMPA)*top_down; 
	duPY_ext_1[thetaa][ii] = -(DT/cmPY)*gGABA*(uPY1[thetaa][ii]-u_GABA)*delta_T*rEXT[thetaa][ii]; 
	if(t>=onset_0 && t<onset_0+period_0 && thetaa<N_assm){    
		I_MGB1[thetaa]  = int_inp0_0*exp(-(pow((thetaa-theta_inp0)/tau_T0,2.0))); 
		I_MGB1[thetaa] += int_inp0_1*exp(-(pow((thetaa-theta_inp1)/tau_T0,2.0)));
		I_MGB1[thetaa] += int_inp0_2*exp(-(pow((thetaa-theta_inp2)/tau_T0,2.0))); 
		I_MGB1[thetaa] += int_inp0_3*exp(-(pow((thetaa-theta_inp3)/tau_T0,2.0)));
		I_MGB1[thetaa] += int_inp0_4*exp(-(pow((thetaa-theta_inp4)/tau_T0,2.0))); 
		I_MGB1[thetaa] += int_inp0_5*exp(-(pow((thetaa-theta_inp5)/tau_T0,2.0))); 
		I_MGB1[thetaa] += int_inp0_6*exp(-(pow((thetaa-theta_inp6)/tau_T0,2.0)));
		I_MGB1[thetaa] += int_inp0_7*exp(-(pow((thetaa-theta_inp7)/tau_T0,2.0)));
	}else{
		I_MGB1[thetaa] = 0.0;
	}
	duPY_MGB1[thetaa] = (DT/cmPY)*(I_MGB1[thetaa]); 

	duPY1 = duPY_leak1[thetaa][ii]                   
	      + duPY_rec_1[thetaa][ii] 
		  + duPY_lat_1[thetaa][ii]  
	      + duPY_topdown[thetaa][ii]
		  + duPY_ext_1[thetaa][ii] 
		  + duPY_MGB1[thetaa]*(int)(rand()/32768.0+inp_prob); 
	return(duPY1);
}
double difuPY2(int thetaa, int ii){
	int jj,thetdash;
	double duPY2=0.0,rec_exT2=0.0,lat_exT=0.0,lat_ihTT=0.0,lat_ihTP=0.0,exP=0.0,c0=0.5,c1=1.0,pw=2.0;
	double bottom_up = 0.0;
	double alpha_T=1.0, tau_T=1.0;   

	duPY_leak2[thetaa][ii] = -(DT*gmPY/cmPY)*(uPY2[thetaa][ii]-UPYres); 
	for (jj=0; jj<N_T; ++jj){ 
		rec_exT2 += w_rec_2[thetaa][ii][thetaa][jj]*rPY2[thetaa][jj]; 
	}
	duPY_rec_2[thetaa][ii] = -(DT/cmPY)*gAMPA*(uPY2[thetaa][ii]-u_AMPA)*rec_exT2; 

	for (thetdash=0; thetdash<=N_assm-1; ++thetdash){
		for (jj=0; jj<=N_T; ++jj){
			bottom_up += w_v2_v1[thetaa][ii][thetdash][jj]*rPY1d[thetdash][jj][t-delay];
		}
	}
	duPY_bottomup[thetaa][ii] = -(DT/cmPY)*gAMPA*(uPY2[thetaa][ii]-u_AMPA)*bottom_up; 
	for (jj=0; jj<N_T; ++jj){ 
		lat_ihTT += w_lat_2[thetaa][ii][jj]*rLB2[thetaa][jj];
	}
	duPY_lat_2[thetaa][ii] = -(DT/cmPY)*gGABA*(uPY2[thetaa][ii]-u_GABA)*lat_ihTT; 

	if(t>=onset_V2 && t<onset_V2+period_V2){         
		I_V2[0] = int_inpV2; 
	}else{
		I_V2[0] = 0.0;
	}
	duPY_inp = (DT/cmPY)*(I_V2[0]); 
	duPY2 = duPY_leak2[thetaa][ii]                  
	      + duPY_rec_2[thetaa][ii] 
		  + duPY_fed_2[thetaa][ii] 
		  + duPY_lat_2[thetaa][ii]  
		  + duPY_bottomup[thetaa][ii]
		  + duPY_ext_2
		  + duPY_inp; 
	return(duPY2);
}
double difuSB1(int thetaa, int ii){
	double duSB1; 

	duSB_leak1[thetaa][ii] = -(DT*gmSB/cmSB)*(uSB1[thetaa][ii]-USBres);     
	if (thetaa<N_assm){
		duSB_1[thetaa][ii] = -(DT/cmSB)*gAMPA*(uSB1[thetaa][ii]-u_AMPA)
						   *(wSB_PY1[thetaa][ii]*rPY1[thetaa][ii]); 
	}
	if (thetaa==N_assm){
		duSB_1[thetaa][ii] = -(DT/cmSB)*gAMPA*(uSB1[thetaa][ii]-u_AMPA)
						   *(wSB_PY1[thetaa][ii]*rPY1[thetaa][ii]); 
	}

	duSB1 = duSB_leak1[thetaa][ii]                   
		  + duSB_1[thetaa][ii]
	      + duSB_ext_1;
	return(duSB1);
}
double difuLB1(int thetaa, int ii){
	int jj,thetdash;
	double duLBB,lat_ih=0.0;
	
	duLB_leak1[thetaa][ii] = -(DT*gmLB/cmLB)*(uLB1[thetaa][ii]-ULBres); 
	for (thetdash=0; thetdash<=N_assm-1; ++thetdash){  
		for (jj=0; jj<=N_T; ++jj){
			lat_ih += (wLB_PY1[thetaa][thetdash][ii]*rPY1[thetdash][jj] + wLB_PY1[thetaa][N_assm][ii]*rPY1[N_assm][jj]); 
		}
	}
	duLB1[thetaa][ii] = -(DT/cmLB)*gAMPA*(uLB1[thetaa][ii]-u_AMPA)*lat_ih; 
	duLBB = duLB_leak1[thetaa][ii]                   
		  + duLB1[thetaa][ii];   
	return(duLBB);
}
double difuLB2(int thetaa, int ii){
	int jj,thetdash;
	double duLBB,lat_ih=0.0;
	
	duLB_leak2[thetaa][ii] = -(DT*gmLB/cmLB)*(uLB2[thetaa][ii]-ULBres); 
	for (thetdash=0; thetdash<=N_assm-1; ++thetdash){  
		for (jj=0; jj<=N_T; ++jj){
			lat_ih += (wLB_PY2[thetaa][thetdash][ii]*rPY2[thetdash][jj] + wLB_PY2[thetaa][N_assm][ii]*rPY2[N_assm][jj]); 
		}
	}
	duLB2[thetaa][ii] = -(DT/cmLB)*gAMPA*(uLB2[thetaa][ii]-u_AMPA)*lat_ih; 
	duLBB = duLB_leak2[thetaa][ii]                   
		  + duLB2[thetaa][ii]   
		  + duLB_ext2;
	return(duLBB);
}
double difuG(int thetaa, int ii){ 
	double duG=0.0,IatoG=0.0,GtoG=0.0;

	duG_leak[thetaa][ii] = -(DT*gmG/cmG)*(uG[thetaa][ii]-UGres); 
    
	IatoG = wG_Ia[thetaa]*(sF[thetaa][ii]*sF[thetaa][ii]*sF[thetaa][ii]*sF[thetaa][ii])/(sF[thetaa][ii]*sF[thetaa][ii]*sF[thetaa][ii]*sF[thetaa][ii]+Kd); 
	duG_Ia[thetaa][ii] = -(DT/cmG)*gGABAb*(uG[thetaa][ii]-u_GABAb)*IatoG; 

	if (ii == 0)  duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][0]-uG[thetaa][19])  + (uG[thetaa][0]-uG[thetaa][1])); 
	if (ii == 1)  duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][1]-uG[thetaa][0])   + (uG[thetaa][1]-uG[thetaa][2])); 
	if (ii == 2)  duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][2]-uG[thetaa][1])   + (uG[thetaa][2]-uG[thetaa][3])); 
	if (ii == 3)  duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][3]-uG[thetaa][2])   + (uG[thetaa][3]-uG[thetaa][4])); 
	if (ii == 4)  duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][4]-uG[thetaa][3])   + (uG[thetaa][4]-uG[thetaa][5])); 
	if (ii == 5)  duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][5]-uG[thetaa][4])   + (uG[thetaa][5]-uG[thetaa][6])); 
	if (ii == 6)  duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][6]-uG[thetaa][5])   + (uG[thetaa][6]-uG[thetaa][7])); 
	if (ii == 7)  duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][7]-uG[thetaa][6])   + (uG[thetaa][7]-uG[thetaa][8]));
	if (ii == 8)  duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][8]-uG[thetaa][9])   + (uG[thetaa][8]-uG[thetaa][10])); 
	if (ii == 9)  duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][9]-uG[thetaa][8])   + (uG[thetaa][9]-uG[thetaa][10])); 
	if (ii == 10) duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][10]-uG[thetaa][9])  + (uG[thetaa][10]-uG[thetaa][11])); 
	if (ii == 11) duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][11]-uG[thetaa][10]) + (uG[thetaa][11]-uG[thetaa][12])); 
	if (ii == 12) duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][12]-uG[thetaa][11]) + (uG[thetaa][12]-uG[thetaa][13])); 
	if (ii == 13) duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][13]-uG[thetaa][12]) + (uG[thetaa][13]-uG[thetaa][14])); 
	if (ii == 14) duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][14]-uG[thetaa][13]) + (uG[thetaa][14]-uG[thetaa][15])); 
	if (ii == 15) duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][15]-uG[thetaa][14]) + (uG[thetaa][15]-uG[thetaa][16])); 
	if (ii == 16) duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][16]-uG[thetaa][15]) + (uG[thetaa][16]-uG[thetaa][17])); 
	if (ii == 17) duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][17]-uG[thetaa][16]) + (uG[thetaa][17]-uG[thetaa][18])); 
	if (ii == 18) duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][18]-uG[thetaa][17]) + (uG[thetaa][18]-uG[thetaa][19])); 
	if (ii == 19) duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][19]-uG[thetaa][18]) + (uG[thetaa][19]-uG[thetaa][0])); 

	duG = duG_leak[thetaa][ii]   
		+ duG_Ia[thetaa][ii] 
		+ duG_G[thetaa][ii];	  
	return(duG);
}

double difrPY1(int thetaa, int ii){
	double drPYT;
	drPYT = DT*(alph_AMPA*GlutPY_c1[thetaa][ii]*(1.0-rPY1[thetaa][ii])-beta_AMPA*rPY1[thetaa][ii]); 
	return(drPYT);
}

double difrPY2(int thetaa, int ii){
	double drPYT;
	drPYT = DT*(alph_AMPA*GlutPY_c2[thetaa][ii]*(1.0-rPY2[thetaa][ii])-beta_AMPA*rPY2[thetaa][ii]); 
	return(drPYT);
}

double difrSB1(int thetaa, int ii){
	double drSBT;
	drSBT = DT*(K1*GABASB_c1[thetaa][ii][t]*(1.0-rSB1[thetaa][ii])-K2*rSB1[thetaa][ii]); 
	return(drSBT);
}
double difsF(int thetaa, int ii){
	double dsF;
	dsF = DT*(K3*rSB1[thetaa][ii]-K4*sF[thetaa][ii]); 
	return(dsF);
}
double difrLB1(int thetaa, int ii){
	double drLB;
	drLB = DT*(alph_GABA*GABALB_c1[thetaa][ii][t]*(1.0-rLB1[thetaa][ii])-beta_GABA*rLB1[thetaa][ii]); 
	return(drLB);
}

double difrLB2(int thetaa, int ii){
	double drLB;
	drLB = DT*(alph_GABA*GABALB_c2[thetaa][ii][t]*(1.0-rLB2[thetaa][ii])-beta_GABA*rLB2[thetaa][ii]); 
	return(drLB);
}

double difrEXT(int thetaa, int ii){
	double drEXT;
	drEXT = DT*(alph_GABA*GABA_ext[thetaa][ii]*(1.0-rEXT[thetaa][ii])-beta_GABA*rEXT[thetaa][ii]); 
	return(drEXT);
}

void dfsGABA_ext(void){
	int theta,i;

	for (theta=0; theta<=N_assm; ++theta){
			for (i=0; i<=N_T; ++i){	
				I_GABA[theta][i] = 0.0;
			}
	}
	for (theta=0; theta<=N_assm; ++theta){
			for (i=0; i<=N_T; ++i){	
				I_GABA[theta][i] += -gamma_*(GABA_ext[theta][i]-GABA_c0)*DT + m_G*(uG[theta][i]-uG_trn)*(GABAamb_max-GABA_ext[theta][i])*(GABA_ext[theta][i]-GABAamb_min)*DT; 
				GABA_extw[theta][i] = GABA_ext[theta][i] + I_GABA[theta][i];
			}
	}
	for (theta=0; theta<=N_assm; ++theta){
			for (i=0; i<=N_T; ++i){	
				GABA_ext[theta][i] = GABA_extw[theta][i];
			}
	}
}

double sigmoidPY(double uPYY){
	return(1/(1+exp(-steep_PY*(uPYY-thres_PY))));
}

double sigmoidPY2(double uPYY){
	return(1/(1+exp(-steep_PY2*(uPYY-thres_PY2))));
}

double sigmoidSB(double uSBB){
	return(1/(1+exp(-steep_SB*(uSBB-thres_SB))));
}

double sigmoidLB(double uLBB){
	return(1/(1+exp(-steep_LB*(uLBB-thres_LB))));
}

double sigmoidLB2(double uLBB){
	return(1/(1+exp(-steep_LB2*(uLBB-thres_LB2))));
}

double rand01(long int *ix){
   double x;
   long int of;
   *ix=(*ix)*48828125;
   if(*ix<0){
	   of=2147483647;
       *ix=(*ix+of)+1;
   }
   x=(double)(*ix)*0.4656613e-9;
   return(x);
}

void init(void){
	int theta,i,thetaa,ii,thetdash,jj;
    double alph_w_L_TV1=2.0; // P-to-BL weight in Nsen    
    double alph_w_L_TV2=2.0; // Nmot 
	double w_v2_1=5.0;       // P(Nsen)-toP(Nmot) weight
	double w_v1P_2=5.0;      // P(Nmot)-toP(Nsen) weight

	srand(17);
	for (thetaa=0; thetaa<=N_assm; ++thetaa){
			for (ii=0; ii<=N_T; ++ii){
     			for (thetdash=0; thetdash<=N_assm; ++thetdash){
					for (jj=0; jj<=N_T; ++jj){
						w_rec_1[N_assm][ii][thetdash][jj] = 0.0; // 
						w_rec_2[N_assm][ii][thetdash][jj] = 0.0; // 
					}
				}
			}
	}

	for (thetaa=0; thetaa<=N_assm; ++thetaa){
			for (ii=0; ii<=N_T; ++ii){
     			for (thetdash=0; thetdash<=N_assm; ++thetdash){
					for (jj=0; jj<=N_T; ++jj){
						w_v2_v1[thetaa][ii][thetdash][jj] = 0.0; // 
					}
				}
			}
	}

	for (thetaa=0; thetaa<=N_assm-1; ++thetaa){  
		for (ii=0; ii<=N_T; ++ii){
			for (jj=0; jj<=N_T; ++jj){
				if (ii!=jj)	w_rec_1[thetaa][ii][thetaa][jj]  = wLPPV1;   
				if (ii!=jj)	w_rec_2[thetaa][ii][thetaa][jj]  = wLPPV2; 
			}
		}
	}

	for (thetaa=0; thetaa<=N_assm-1; ++thetaa){  
		for (ii=0; ii<=N_T; ++ii){
			for (jj=0; jj<=N_T; ++jj){
				w_lat_1[thetaa][ii][jj] = 4.0;  
				w_lat_2[thetaa][ii][jj] = 4.0;     
			}
		}
	}

	for (thetaa=0; thetaa<=N_assm-1; ++thetaa){  
		for (ii=0; ii<=N_T; ++ii){
			wSB_PY1[thetaa][ii] = 60.0;  
		}
	}

	for (thetaa=0; thetaa<=N_assm; ++thetaa){
		for (thetdash=0; thetdash<=N_assm; ++thetdash){
			for (ii=0; ii<=N_T; ++ii){
				if (thetaa != thetdash){
					wLB_PY1[thetaa][thetdash][ii] = alph_w_L_TV1; 
					wLB_PY2[thetaa][thetdash][ii] = alph_w_L_TV2; 
				}
			}
		}
	}

	for (thetaa=0; thetaa<=N_assm; ++thetaa){  
		wG_Ia[thetaa] = 3.0; 
	}

	for (thetaa=0; thetaa<=N_assm-1; ++thetaa){  
		for (ii=0; ii<=N_T; ++ii){
			for (jj=0; jj<=N_T; ++jj) w_v1P_v2[thetaa][ii][thetaa][jj] = w_v1P_2;
		}
	}

	for (thetaa=0; thetaa<=N_assm-1; ++thetaa){  
		for (ii=0; ii<=N_T; ++ii){
			for (jj=0; jj<=N_T; ++jj) w_v2_v1[thetaa][ii][thetaa][jj] = w_v2_1;
		}
	}
	for (theta=0; theta<=N_assm; ++theta){
		for (i=0; i<=N_T; ++i){	
			GABA_ext[theta][i]  = GABA_c0; 
			GABA_extw[theta][i] = GABA_c0; 
		}
	}
}
