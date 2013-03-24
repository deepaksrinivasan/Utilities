#include "mex.h"
#include <matrix.h>
#include <math.h>




double compat(int a,int i,int b,int j,double g1,double g2)
{
double alpha=0.1; double c;
if (((a-i)*(b-j)) < 0)
                   {c=0;}
else
                    
                    c=g1*g2/(1+alpha*abs(abs(a-b)-abs(i-j)));
                    //c=(1-3*abs(g1-g2));%/max(abs(a-b)-abs(i-j),1); %penalizes short range long range edge mapping

return c;

}

void eye_init(int num_nodes1,int num_nodes2,double epsilon,double *M)
{
	int i,j;
	for (i=0;i<num_nodes1;i++)
	{
		for (j=0;j<num_nodes2;j++)
			{
				if (i==j)
					M[(i*num_nodes2)+j]=1+epsilon;
				else
					M[(i*num_nodes2)+j]=1;
				//mexPrintf("%f\t",M[i*num_nodes2+j]);
		
			}
			//mexPrintf("\n");
	}
}

void fill_arrayformatlab(int num_nodes1, int num_nodes2, double *M1,double *M)
{
	int i, j;
	for(i=0;i<num_nodes1;i++)
		{
		for(j=0;j<num_nodes2;j++)
			{
			M1[j*num_nodes1+i] = M [i*num_nodes2+j];
			}
		}
}

void fill_arrayforc(int num_nodes1, int num_nodes2, double *M1,double *M)
{
	int i, j;
	for(i=0;i<num_nodes1;i++)
		{
		for(j=0;j<num_nodes2;j++)
			{
			M1[i*num_nodes2+j] = M [j*num_nodes1+i];
			}
		}
}


			











void mexFunction(int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
int m, n;
double *G, *g, *M, *Q, *M1ptr, *G_ones, *g_ones;
double *d1, *d2 ; 
double c;
double beta_i=0.5, beta_f=10, neta=1.075,epsilon=0.005, n_iter=10;
int n1,a,i,b,j,p,p1,p2, b1, j1,i2,j2;
mxArray *G_sum[1], *g_sum[1], *M1[1], *M2[1];
int num_nodes1=mxGetNumberOfElements(prhs[0]);
int num_nodes2=mxGetNumberOfElements(prhs[1]);
mxArray *nodes1[1], *nodes2[1];


plhs[0]=mxCreateDoubleMatrix(num_nodes1,num_nodes2,mxREAL);
plhs[1]=mxCreateDoubleMatrix(num_nodes1,num_nodes2,mxREAL);
M=mxGetPr(plhs[0]);
Q=mxGetPr(plhs[1]);
G=mxGetPr(prhs[2]);
g=mxGetPr(prhs[3]);


 M1[0] = mxCreateDoubleMatrix(num_nodes1, num_nodes2,mxREAL);
M2[0] = mxCreateDoubleMatrix(num_nodes1, num_nodes2,mxREAL);
M1ptr=mxGetPr(M1[0]);
 

eye_init(num_nodes1,num_nodes2,epsilon,M);
fill_arrayformatlab(num_nodes1,num_nodes2,mxGetPr(M1[0]),M);
fill_arrayformatlab(num_nodes1,num_nodes2,mxGetPr(M2[0]),M);
fill_arrayforc(num_nodes1,num_nodes2,M,mxGetPr(M2[0]));
//fill_arrayformatlab(num_nodes1,num_nodes1,G,G);
//fill_arrayformatlab(num_nodes2,num_nodes2,g,g);

/*mexCallMATLAB(0,NULL,1,M1,"disp");
mexCallMATLAB(0,NULL,1,M2,"disp");






for (i=0;i<num_nodes1;i++)
	{
	for (j=0;j<num_nodes2;j++)
		{
		
		mexPrintf("%f\t",M[(i*num_nodes2)+j]);
		}
	mexPrintf("\n");
	}
*/


double beta=beta_i;int dummy=1;
while (beta<=beta_f)
	{
	for(n1=0;n1<4;n1++)
		{
		for (a=0;a<num_nodes1;a++)
			{
			for (i=0;i<num_nodes2;i++)
				{
				Q[a*num_nodes2+i]=0;
				
				nodes1[0]=mxGetFieldByNumber(prhs[0],a,0);
                                
	
		
				d1=mxGetPr(nodes1[0]);
                            

				for (b=0;b<mxGetNumberOfElements(mxGetFieldByNumber (prhs[0], a, 0));b++)
					{
					nodes2[0]=mxGetFieldByNumber(prhs[1],i,0);
					d2=mxGetPr(nodes2[0]);
                                        
                                        

					for (j=0;j<mxGetNumberOfElements(mxGetFieldByNumber (prhs[1], i, 0));j++)
						
						{
	
                                             
						
						
						c=compat(a,i,b,j,G[(a*num_nodes1)+ ((int) d1[b]-1)],g[(i*num_nodes2)+ ((int) d2[j]-1)]);
						//c=compat(a,i,b,j,1,1);
						Q[a*num_nodes2+i]+=(M[b*num_nodes2+j]*c);
						
						
						}
					
					}
				
				
				
				}
			}
	
                
	
			
		for ( p1=0;p1<num_nodes1;p1++)
			{
			for ( p2=0;p2<num_nodes2;p2++)
				{
				M[p1*num_nodes2+p2]=exp(beta*Q[p1*num_nodes2+p2]);
					
				}
				
			}
			
			
		//memcpy(mxGetPr(M1[0]),M, num_nodes1*num_nodes2*sizeof(double));
		fill_arrayformatlab(num_nodes1,num_nodes2,mxGetPr(M1[0]),M);
		
		mexCallMATLAB(1,M2,1,M1,"rowcol_normal");
		 
		
		
		fill_arrayforc(num_nodes1,num_nodes2,M,mxGetPr(M2[0]));

		
		
		
		
		}
	beta*=neta;
	//mexPrintf("%f\n",beta);
	//mexCallMATLAB(0,NULL,1,nodes1,"disp");
	
	//mexCallMATLAB(0,NULL,1,nodes2,"disp");
	
	
	}

/*for (a=0;a<num_nodes1;a++)
			{for (b=0;b<mxGetNumberOfElements(mxGetFieldByNumber (prhs[0], a, 0));b++)
                           {
	
                                 nodes1[0]=mxGetFieldByNumber(prhs[0],a,0);
                                
	
		
				d1=mxGetPr(nodes1[0]);
                                mexPrintf("%f\t",G[(a*num_nodes1)+ ((int) d1[b]-1)]);

}
                                mexPrintf("\n");
}
mexPrintf("\n");
for (i=0;i<num_nodes2;i++)
			{for (j=0;j<mxGetNumberOfElements(mxGetFieldByNumber (prhs[1], i, 0));j++)
                           {
	
                                 nodes2[0]=mxGetFieldByNumber(prhs[1],i,0);
                                
	
		
				d2=mxGetPr(nodes2[0]);
                                mexPrintf("%f\t",g[(i*num_nodes2)+ ((int) d2[j]-1)]);

}
                                mexPrintf("\n");
}*/
}





	
			
		
		
		


