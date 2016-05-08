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
	std::clock_t start;
   	double duration;
	std::cout<<"\n\n\n\n\n\n\nStart again\n\n";
    	start = std::clock();
    	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	start = std::clock();	
	//Preprocessing colours in vector col
	/*{{{*/
	double4_t* sRect = double4_alloc((nx+1)*(ny+1));
	double4_t totBack;
	double4_t* col = double4_alloc(((nx+1)*(ny+1)));
	travx = nx+1;
	totBack = double4_0;
	std::cout<<"Colours:  "<<std::endl;
	for(i = 0; i<= ny ; ++i){
		for(j = 0; j<= nx; ++j){
			//std::cout<<"check this one ("<<i<<","<<j<<")"<<"  "<<std::endl;
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
			std::cout<<"("<<j<<","<<i<<")"<<col[i*travx+j][0]<<"  ";
		}	
		std::cout<<std::endl;
	}
	std::cout<<std::endl<<"col again"<<std::endl;
;
	for(i =0; i<=ny;i++){
		for(j=0;j<=nx;j++){
	std::cout<<j<<","<<i<<"    "<<col[i*travx+j][0]<<"  ";
		}	
		std::cout<<std::endl;
	}
	
	std::cout<<std::endl<<"SUM CALCULATION sRect additions"<<std::endl;
	//making sum from 0,0 point to x,y and saving them in sReact
	/*{{{*/

	for( i = 1; i<= ny; ++i){
		for( j= 1; j<= nx; ++j){
			std::cout<<std::endl<<"("<<j<<","<<i<<")"<<std::endl;
			sRect[i*travx+j] = double4_0;
			for( ii = i; ii>=0 ; --ii){
				for(jj = j; jj>=0; --jj){
					sRect[i*travx+j] += col[ii*travx + jj];
					std::cout<<"("<<jj<<","<<ii<<") "<<col[ii*travx + jj][0]<<" + ";
				}				
			}
			std::cout<<std::endl;
		}			
	}

	std::cout<<std::endl<<"Check for sRect pos and values"<<std::endl;
	std::cout<<sRect[3*travx+1][0]<<"   tot is  "<<totBack[0]<<std::endl;
	for( i = 0 ; i <=ny ; i++){
		for(j=0 ; j <=nx ; j++){
			std::cout<<"("<<j<<","<<i<<")"<<sRect[i*travx + j][0]<<"  ";
		}
		std::cout<<std::endl;
	}
	/*}}}*/	
	double maxVal = 0;
	double tempMax = 0;
	double X;		//reciprocal of X
	double Y;		//reciprocal of Y
	int totSize = nx*ny;
//	int K,M;
	double4_t Hxy = double4_0;
	//traversing through pixels with diff size of boxes
	//making the box 2 consecutive for loops (0,0) constant
	// i*j is the size of the boxes
	for(i = 1; i<=ny; ++i){
		for(j=1; j<=nx; ++j){
				
				std::cout<<std::endl<<" i is "<<i<<" j is "<<j<<std::endl;
				std::cout<<"--------------------------------------"<<std::endl;
				ii = i*j;
				X = 1.0/double(ii);
				jj = totSize - ii;
				Y = 1.0/double(jj);
				
				
				//moving the box 2 consecutive for loops
				for( k= 0; k<ny-i+1; ++k){
					 double4_t Vxy = double4_0;
					 double4_t Yxy = double4_0;
					for(m = 0;m<nx-j+1; ++m){
							std::cout<<"value of x y are "<<m<<"  "<<k<<std::endl;
							std::cout<<"values inside rect are "<<m+j<<"  "<<k+i<<std::endl;
							//std::cout<<"value to check is "<<sRect[(k+i)*travx + (m+j)][0]<<std::endl;
							//std::cout<<"value to check is "<<sRect[(k)*travx + (m+j)][0]<<std::endl;
							if( k ==0 && m==0){
								Vxy = sRect[(k+i)*travx + (m+j)];		
							}
							else if( m==0){
								Vxy =sRect[(k+i)*travx + (m+j)]-sRect[(k)*travx + (m+j)];
							}
							else if( k==0){
								Vxy = sRect[(k+i)*travx + (m+j)]-sRect[(k+i)*travx + (m)];
							}	
							else{
								Vxy = sRect[(k+i)*travx + (m+j)] -sRect[(k+i)*travx + (m)]-sRect[(k)*travx + (m+j)]+sRect[(k)*travx + (m)];
							}
							Yxy = totBack - Vxy;
							Hxy[0] = pow(Vxy[0],2)*X + pow(Yxy[0],2)*Y;
							Hxy[1] = pow(Vxy[1],2)*X + pow(Yxy[1],2)*Y;
							Hxy[0] = pow(Vxy[2],2)*X + pow(Yxy[2],2)*Y;		
							tempMax = Hxy[0] + Hxy[1] + Hxy[2];
							if(tempMax >= maxVal){
								std::cout<<std::endl<<"y0 x0 y1 x1 are "<<k<<"  "<<m<<"  "<<k+i<<"  "<<m+j<<std::endl;
								maxVal = tempMax;
								result.x0 = m;
								result.y0 = k;
								result.x1 = m+j;
								result.y1 = k+i;
								result.inner[0] = col[(k+i)*travx + (m+j)][0];
								result.inner[1] = col[(k+i)*travx + (m+j)][1];
								result.inner[2] = col[(k+i)*travx + (m+j)][2];
							}				
//							std::cout<<std::endl<<"y0 x0 y1 x1 are "<<result.y0<<"  "<<result.x0<<"  "<<result.y1<<"  "<<result.x1<<std::endl;
//							std::cout<<std::endl<<"colors of inner are  "<<result.inner[0]<<result.inner[1]<<result.inner[2]<<std::endl;
						}
					}
					std::cout<<std::endl;
				}
		}
	
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        std::cout<<"time spent for whole : "<<(float) duration <<'\n';
        std::cout<<std::endl;

	return result;
}
