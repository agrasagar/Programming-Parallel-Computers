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
	travx = (nx + 1);
	totBack = double4_0;
	std::cout<<std::endl;
	for(i = 0; i< ny+1 ; ++i){
		for(j = 0; j< nx+1; ++j){
			col[i*travx+j] = double4_0;
			if( j!=0){
				for(c = 0; c<3; ++c){
					col[i*travx+j][c] = data[c + 3 * (j-1) + 3 * nx * (i)];
				}
				totBack+=col[i*travx+j];
			}
			std::cout<<"("<<i<<","<<j<<")"<<col[i*travx+j][0]<<"  ";
		}	
		std::cout<<std::endl;
	}
//	std::cout<<std::endl<<"SUM CALCULATION"<<std::endl;
	//making sum from 0,0 point to x,y and saving them in sRectc
	/*{{{*/
	for( i = 0; i< ny+1; ++i){
		for( j= 1; j< nx+1; ++j){
//			std::cout<<std::endl<<"("<<i<<","<<j<<")"<<std::endl;
			sRect[i*travx+j] = double4_0;
			for( ii = i; ii>=0 ; --ii){
				for(jj = j; jj>=0; --jj){
					sRect[i*travx+j] += col[ii*travx + jj];
//					std::cout<<std::endl<<"("<<ii<<","<<jj<<") "<<col[ii*travx + jj][0]<<" + ";
				}				
			}
//			std::cout<<std::endl;
		}			
	}

//	std::cout<<std::endl<<"Check for total"<<std::endl;
//	std::cout<<sRect[2*travx+1][0]<<"   tot is  "<<totBack[0]<<std::endl;
	for( i = 0 ; i < ny+1 ; i++){
		for(j=1 ; j < nx+1 ; j++){
//			std::cout<<std::endl<<"value of srect 0 "<<sRect[i*travx + j][0]<<std::endl;
		}
	}
	/*}}}*/	
	double maxVal = 0;
	double tempMax = 0;
	double X;		//reciprocal of X
	double Y;		//reciprocal of Y
	int totSize = nx*ny;
	int K,M;
	double4_t Hxy = double4_0;
	//traversing through pixels with diff size of boxes
	//making the box 2 consecutive for loops (0,0) constant
	// i*j is the size of the boxes
	for(i = 1; i<ny+1; ++i){
		for(j=1; j<nx+1; ++j){
				
				std::cout<<std::endl<<" i is "<<i<<" j is "<<j<<std::endl;
				std::cout<<"--------------------------------------"<<std::endl;
				ii = i*j;
				//if( ii == totSize){
				//	break;
				//}
				X = 1.0/double(ii);
				jj = totSize - ii;
				Y = 1.0/double(jj);
				(void)X;
				(void)Y;
				//moving the box 2 consecutive for loops
				for( k= 0; k<=ny-i; ++k){
					 double4_t Vxy = double4_0;
					 double4_t Yxy = double4_0;
					for(m = 0;m<=nx-j; ++m){
						K = k+i-1;
						M = m+j;
						std::cout<<std::endl<<"value of start x y "<<k<<m<<std::endl;
						std::cout<<std::endl<<"value of end x y "<<K<<M<<std::endl;
						if( m==0 && k==0 ){
							Vxy += sRect[K*travx + M];
							std::cout<<std::endl<<"value of end In m&k==0 loop  x y "<<K<<M<<std::endl;
						}
						else if(m==0){
							Vxy += sRect[K*travx + M]- sRect[(k-1)*travx + (M)];
							std::cout<<std::endl<<"value of end In m==0 loop  x y "<<K<<M<<std::endl;
							std::cout<<"value of end In  loop  K-1 M is  "<<k-1<<M<<std::endl;
						}	
						else if(k==0){
							Vxy += sRect[K*travx + M]- sRect[(K)*travx + (m)];
							std::cout<<std::endl<<"value of end In k==0 loop  x y "<<K<<M<<std::endl;
                                                        std::cout<<"value of end In  loop  K m is  "<<K<<m<<std::endl;
						}			
						else{
							Vxy += sRect[K*travx + M]- sRect[(K)*travx + (m)]- sRect[(k-1)*travx + (M)] + sRect[(k-1)*travx + (m)];
						}
		
						Yxy += totBack-Vxy;
			
						Hxy[0] = pow(Vxy[0],2)*X + pow(Yxy[0],2)*Y;
						Hxy[1] = pow(Vxy[1],2)*X + pow(Yxy[1],2)*Y;
						Hxy[2] = pow(Vxy[1],2)*X + pow(Yxy[2],2)*Y;
						tempMax = Hxy[0]+Hxy[1]+Hxy[2];
						if(tempMax >= maxVal){
							maxVal = tempMax;
							std::cout<<std::endl<<"///////////////////////////////  x0 y0 x1 y1 are "<<k<<"  "<<m<<"  "<<M<<"  "<<K+1<<std::endl;
							result.x0 = m;
							result.y0 = k;
							result.y1 = k+i-1;
							result.x1 = M;
							result.inner[0] = col[K*travx + M][0];
							result.inner[1] = col[K*travx + M][1];
							result.inner[2] = col[K*travx + M][2];
							double val = 1.0/(jj);	
							result.outer[0] = Yxy[0]*val;
							result.outer[1] = Yxy[1]*val;
							result.outer[2] = Yxy[2]*val;
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
