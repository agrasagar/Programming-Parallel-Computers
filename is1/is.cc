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
	
	//sum calculation
	for( i = 1; i<= ny; ++i){
		for( j = 1; j<= nx ; ++j){
			sRect[i*travx + j] += col[(i)*travx + (j)] + sRect[(i-1)*travx + (j)] + sRect[(i)*travx + (j-1)] - sRect[(i-1)*travx + (j-1)];
		}
	}
	

	
	double maxVal = 0;
	int totSize = nx*ny;
	
	{
	//traversing through pixels with diff size of boxes
	//making the box 2 consecutive for loops (0,0) constant
	#pragma omp parallel for schedule(dynamic) private(k,m,i,j)
	for(i = 1; i<=ny; ++i){
		for(j=1; j<=nx; ++j){
				
				int totS = i*j;
				double X = 1.0/double(totS);
				int  totR = totSize - totS;
				double Y = 1.0/double(totR);
				int condk = ny-i+1;
				int condm = nx-j+1;	
				double tempMax = 0;
				double maxValPriv = 0;
				int tmpm = 0,tmpk = 0;
                                double4_t tmpVxy = double4_0;
                                double4_t tmpYxy = double4_0;
				//moving the box 2 consecutive for loops
				for( k= 0; k<condk; ++k){
					double4_t Hxy = double4_0;
					double4_t Vxy = double4_0;
                                	double4_t Yxy = double4_0;
					for(m = 0;m<condm; ++m){
							
						Vxy = sRect[(k+i)*travx + (m+j)] -sRect[(k+i)*travx + (m)]-sRect[(k)*travx + (m+j)]+sRect[(k)*travx + (m)];
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
							}				
						}
					}
					#pragma omp critical
					if(maxValPriv > maxVal){
                                        	maxVal = maxValPriv;
                                                result.x0 = tmpm;
                                                result.y0 = tmpk;
                                                result.x1 = tmpm+j;
                                                result.y1 = tmpk+i;
                                                result.inner[0] = tmpVxy[0]*X;
                                                result.inner[1] = tmpVxy[1]*X;
                                                result.inner[2] = tmpVxy[2]*X;
                                                result.outer[0] = tmpYxy[0]*Y;
                                                result.outer[1] = tmpYxy[1]*Y;
                                                result.outer[2] = tmpYxy[2]*Y;
                        
                                       }      					

			}
		}
	}
	
	return result;
}
