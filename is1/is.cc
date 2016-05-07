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
    	/* Your algorithm here */
    	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    	//std::cout<<"printf: "<< duration <<'\n';
	start = std::clock();	
	//Preprocessing colours in vector col
	/*{{{*/
	double4_t* sRect = double4_alloc(nx*ny);
	double4_t totBack;
	double4_t* col = double4_alloc(nx*ny);
	travx = 0;
	totBack = double4_0;

	for(i = 0; i<ny ; ++i){
		for(j = 0; j<nx; ++j){
			col[travx] = double4_0;
			for(c = 0; c<3; ++c){
				col[travx][c] = data[c + 3 * j + 3 * nx * i];
			}
			totBack+=col[travx];
			travx++;
		}	
	}
	//making sum from 0,0 point to x,y and saving them in sRectc
	/*{{{*/
	travx = nx;
	for( i = 0; i<ny; ++i){
		for( j= 0; j< nx; ++j){
			sRect[i*travx+j] = double4_0;
			for( ii = i; ii>=0 ; --ii){
				for(jj = j; jj>=0; --jj){
					sRect[i*travx+j] += col[ii*travx + jj];
				}				
			}
		}			
	}
	/*}}}*/	
	double maxVal = 0;
	double tempMax = 0;
	int X;		//reciprocal of X
	int Y;		//reciprocal of Y
	int totSize = nx*ny;
	double4_t Hxy = double4_0;
	//traversing through pixels with diff size of boxes
	//making the box 2 consecutive for loops (0,0) constant
	// i*j is the size of the boxes
	for(i = 1; i<=ny; ++i){
		for( j=1; j<=nx; ++j){
				X = 1/(i*j);
				Y = 1/(totSize - (i*j)); 		
				//moving the box 2 consecutive for loops
				for( k= 0; k<ny-k+1; ++k){
					 double4_t Vxy = double4_0;
					for(m = 0;m<nx-m+1; ++m){
						if(k == 0 && m == 0){
							Vxy = sRect[(k + i -1)*travx + (m + j - 1)];
						}
						else if(k == 0){
							Vxy = sRect[(k + i -1)*travx + (m + j - 1)] - sRect[(k + i -1)*travx +(m - 1)];
						}
						else if(m == 0){
							Vxy = sRect[(k + i -1)*travx + (m + j - 1)] - sRect[(k - 1)*travx + (m + j - 1)];
						}
						else{
							Vxy = sRect[(k + i -1)*travx + (m + j - 1)] - sRect[(k + i - 1)*travx + (m + j - 2)] - sRect[(k - 1)*travx + (m + j - 1)] + sRect[(k - 1)*travx + (m -1)];
						}
						Hxy[0] = pow(Vxy[0],2)*X + pow(totBack[0],2)*Y;
						Hxy[1] = pow(Vxy[1],2)*X + pow(totBack[1],2)*Y;
						Hxy[2] = pow(Vxy[2],2)*X + pow(totBack[2],2)*Y;
						tempMax = Hxy[0] + Hxy[1] + Hxy[2];
						std::cout<<std::endl<<"Goes Once"<<std::endl;
						std::cout<<std::endl<<"Temp Val = "<<tempMax<<"Max val is= "<<maxVal<<std::endl;
						std::cout<<std::endl<<"-------------------------"<<std::endl;
						if(tempMax > maxVal){
							std::cout<<std::endl<<"Goes Once"<<std::endl;
							maxVal = tempMax;
							result.x0 = (m);
							result.y0 = (k);
							result.x1 = (m + j);
							result.y1 = (k + i);
							result.inner[0] = col[(k + i -1)*travx + (m + j - 1)][0];
							result.inner[1] = col[(k + i -1)*travx + (m + j - 1)][1];
							result.inner[2] = col[(k + i -1)*travx + (m + j - 1)][2];
							if( m-1 >= 0){
								result.outer[0] = col[(k)*travx + (m - 1)][0];
								result.outer[1] = col[(k)*travx + (m - 1)][1];
								result.outer[2] = col[(k)*travx + (m - 1)][2];
							}
							else if( (m+j) < nx ){
								result.outer[0] = col[(k)*travx + (m + j)][0];
								result.outer[1] = col[(k)*travx + (m + j)][1];
								result.outer[2] = col[(k)*travx + (m + j)][2];
							}
							else if( k-1 >= 0){
								result.outer[0] = col[(k-1)*travx + (m)][0];
								result.outer[1] = col[(k-1)*travx + (m)][1];
								result.outer[2] = col[(k-1)*travx + (m)][2];
							}
							else if( (k+i) < ny){	
								result.outer[0] = col[(k+i)*travx + (m)][0];
								result.outer[1] = col[(k+i)*travx + (m)][1];
								result.outer[2] = col[(k+i)*travx + (m)][2];
							}



							std::cout<<std::endl<<"y0 x0 y1 x1 are "<<result.y0<<"  "<<result.x0<<"  "<<result.y1<<"  "<<result.x1<<std::endl;
							std::cout<<std::endl<<"colors are "<<result.inner[0]<<result.inner[1]<<result.inner[2]<<std::endl;
						}
					}
				}
		}
	}
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        std::cout<<"time spent for whole : "<<(float) duration <<'\n';
        std::cout<<std::endl;




	return result;
}
