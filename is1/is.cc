#include "is.h"
#include <stdlib.h>
#include "../common/vector.h"
#include <iostream>
#include <ctime>
#include <math.h>

Result segment(int ny, int nx, const float* data) {
	// FIXME
	Result result { ny/3, nx/3, 2*ny/3, 2*nx/3, {0.0f, 0.0f, 1.0f}, {1.0f, 0.0f, 0.0f} };
	int i,j,ii,jj,k,m,c, travx;
	//Preprocessing colours in vector col
	/*{{{*/
	double4_t* sRect = double4_alloc((nx+1)*(ny+1));
	double4_t totBack;
	double4_t* col = double4_alloc(((nx+1)*(ny+1)));
	travx = nx+1;
	totBack = double4_0;
	for(i = 0; i<= ny ; ++i){
		for(j = 0; j<= nx; ++j){
			if(i == 0 || j ==0){
				 col[i*travx+j] = double4_0;
				 sRect[i*travx+j] = double4_0;
			}
			else{
				col[i*travx+j] = double4_0;
				for(c = 0; c<3; ++c){
					col[i*travx+j][c] = data[c + 3 * (j-1) + 3 * nx * (i-1)];
				}
				totBack+=col[i*travx+j];
			}
			sRect[i*travx+j] = double4_0;
		}	
	}

	for(i =0; i<=ny;i++){
		for(j=0;j<=nx;j++){
		}	
	}
	
	//making sum from 0,0 point to x,y and saving them in sReact
	/*{{{*/

	for( i = 1; i<= ny; ++i){
		for( j= 1; j<= nx; ++j){
			sRect[i*travx+j] = double4_0;
			for( ii = i; ii>=0 ; --ii){
				for(jj = j; jj>=0; --jj){
					sRect[i*travx+j] += col[ii*travx + jj];
				}				
			}
		}			
	}

	 #pragma omp parallel
	{ 
	/*}}}*/	
	double maxVal = 0;
	double tempMax = 0;
	double X;		//reciprocal of X
	double Y;		//reciprocal of Y
	int totSize = nx*ny;
	double4_t Hxy = double4_0;
	//traversing through pixels with diff size of boxes
	//making the box 2 consecutive for loops (0,0) constant
	 #pragma omp for schedule(static,1)
	for(i = 1; i<=ny; ++i){
		for(j=1; j<=nx; ++j){
				
				ii = i*j;
				X = 1.0/double(ii);
				jj = totSize - ii;
				Y = 1.0/double(jj);
				
				
				//moving the box 2 consecutive for loops
				for( k= 0; k<ny-i+1; ++k){
					 double4_t Vxy = double4_0;
					 double4_t Yxy = double4_0;
					for(m = 0;m<nx-j+1; ++m){
//							if( k ==0 && m==0){
//								Vxy = sRect[(k+i)*travx + (m+j)];		
//							}
//							else if( m==0){
//								Vxy =sRect[(k+i)*travx + (m+j)]-sRect[(k)*travx + (m+j)];
//							}
//							else if( k==0){
//								Vxy = sRect[(k+i)*travx + (m+j)]-sRect[(k+i)*travx + (m)];
//							}	
//							else{
							Vxy = sRect[(k+i)*travx + (m+j)] -sRect[(k+i)*travx + (m)]-sRect[(k)*travx + (m+j)]+sRect[(k)*travx + (m)];
//							}
							Yxy = totBack - Vxy;
							Hxy[0] = pow(Vxy[0],2)*X + pow(Yxy[0],2)*Y;
							Hxy[1] = pow(Vxy[1],2)*X + pow(Yxy[1],2)*Y;
							Hxy[2] = pow(Vxy[2],2)*X + pow(Yxy[2],2)*Y;		
							tempMax = Hxy[0] + Hxy[1] + Hxy[2];
							#pragma omp critical
							if(tempMax >= maxVal){
								maxVal = tempMax;
								result.x0 = m;
								result.y0 = k;
								result.x1 = m+j;
								result.y1 = k+i;
								result.inner[0] = Vxy[0]*X;
								result.inner[1] = Vxy[1]*X;
								result.inner[2] = Vxy[2]*X;
								result.outer[0] = Yxy[0]*Y;
								result.outer[1] = Yxy[1]*Y;
								result.outer[2] = Yxy[2]*Y;

							}				
						}
					}
				}
		}
	}
	
	return result;
}
