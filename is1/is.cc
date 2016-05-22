#include "is.h"
#include <stdlib.h>
#include "../common/vector.h"
#include <iostream>
#include <ctime>
#include <math.h>

Result segment(int ny, int nx, const float* data) {
	// FIXME
	Result result { ny/3, nx/3, 2*ny/3, 2*nx/3, {0.0f, 0.0f, 1.0f}, {1.0f, 0.0f, 0.0f} };
	
	int i,j,k,m,c,travx;
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
	
	//sum calculation
	for( i = 1; i<= ny; ++i){
		for( j = 1; j<= nx ; ++j){
			sRect[i*travx + j] += col[(i)*travx + (j)] + sRect[(i-1)*travx + (j)] + sRect[(i)*travx + (j-1)] - sRect[(i-1)*travx + (j-1)];
		}
	}
	

	
	double maxVal = 0;
	int totSize = nx*ny;
	#pragma omp parallel private(k,m,i,j)
	{
	int tmpm = 0,tmpk = 0, tmpi = 0, tmpj = 0;
        double4_t tmpVxy = double4_0;
        double4_t tmpYxy = double4_0;
	double tmpX = 0, tmpY=0;
	double maxValPriv = 0;
	//traversing through pixels with diff size of boxes
	//making the box 2 consecutive for loops (0,0) constant
	#pragma omp for schedule(static,1) 
	for(i = 1; i<=ny; ++i){
		for(j=1; j<=nx; ++j){
				
				int totS = i*j;
				double X = 1.0/double(totS);
				int  totR = totSize - totS;
				double Y = 1.0/double(totR);
				int condk = ny-i+1;
				int condm = nx-j+1;	
				double tempMax = 0;
				int ki,mj;
				//moving the box 2 consecutive for loops
				for( k= 0; k<condk; ++k){
					double4_t Hxy = double4_0;
					double4_t Vxy = double4_0;
                                	double4_t Yxy = double4_0;
					for(m = 0;m<condm; ++m){
						
						ki = k + i;
						mj = m + j;	
						Vxy = sRect[(ki)*travx + (mj)] -sRect[(ki)*travx + (m)]-sRect[(k)*travx + (mj)]+sRect[(k)*travx + (m)];
						Yxy = totBack - Vxy;
							
						Hxy[0] = pow(Vxy[0],2)*X + pow(Yxy[0],2)*Y;
						Hxy[1] = pow(Vxy[1],2)*X + pow(Yxy[1],2)*Y;
						Hxy[2] = pow(Vxy[2],2)*X + pow(Yxy[2],2)*Y;		
						tempMax = Hxy[0] + Hxy[1] + Hxy[2];
							
							if(tempMax >= maxValPriv){
								maxValPriv = tempMax;
								tmpm = m;
								tmpk = k;	
								tmpVxy = Vxy;	
								tmpYxy = Yxy;
								tmpi = i;
								tmpj = j;
								tmpX = X;
								tmpY = Y;
							}				
						}
					}
			}
		}
		#pragma omp critical
                if(maxValPriv > maxVal){
                	maxVal = maxValPriv;
                        result.x0 = tmpm;
                        result.y0 = tmpk;
                        result.x1 = tmpm+tmpj;
                        result.y1 = tmpk+tmpi;
                        
			result.inner[0] = (float)tmpVxy[0]*tmpX;
                        result.inner[1] = (float)tmpVxy[1]*tmpX;
                        result.inner[2] = (float)tmpVxy[2]*tmpX;
                        result.outer[0] = (float)tmpYxy[0]*tmpY;
                        result.outer[1] = (float)tmpYxy[1]*tmpY;
                        result.outer[2] = (float)tmpYxy[2]*tmpY;
                        
               }                                        
	}
	
	return result;
}
