#include "is.h"
#include <stdlib.h>
#include "../common/vector.h"
#include <iostream>
#include <ctime>

Result segment(int ny, int nx, const float* data) {
	// FIXME
	Result result { ny/3, nx/3, 2*ny/3, 2*nx/3, {0.0f, 0.0f, 1.0f}, {1.0f, 0.0f, 0.0f} };
	int i,j,jj,ii,c, travx;
	std::clock_t start;
   	double duration;

    	start = std::clock();
    	/* Your algorithm here */
    	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    	//std::cout<<"printf: "<< duration <<'\n';
	start = std::clock();	
	//Preprocessing colours in vector col
	/*{{{*/
	double* sRectc1 = (double*)malloc((nx*ny)*sizeof(double));
	double* sRectc2 = (double*)malloc((nx*ny)*sizeof(double));
	double* sRectc3 = (double*)malloc((nx*ny)*sizeof(double));

	double4_t totBack;
	double4_t* col = double4_alloc(nx*ny);
	travx = 0;
	totBack = double4_0;

	for(i = 0; i<ny ; ++i){
		for(j = 0; j<nx; ++j){
			col[travx] = double4_0;
			for(c = 0; c<3; ++c){
				col[travx][c] = data[c + 3 * j + 3 * nx * i];
				std::cout<<std::endl<<"colour value is "<<data[c + 3 * j + 3 * nx * i];
			}
			totBack+=col[travx];
			travx++;
			std::cout<<std::endl;
		}	
	}
	std::cout<<std::endl;
	for(i =0; i< (nx*ny); i++){
		std::cout<<col[i][0]<<" "<<col[i][1]<<" "<<col[i][2]<<" "<<col[i][3];
		std::cout<<std::endl;
	}
	std::cout<<"total is "<<totBack[0]<<" "<<totBack[1]<<" "<<totBack[2]<<" "<<totBack[3]<<std::endl;	
	/**}}}*/
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        std::cout<<"time spent: "<< duration <<'\n';
	std::cout<<std::endl;
	//making sum from 0,0 point to x,y and saving them in sRectc(n)
	/*{{{*/
	for( i = 0; i<ny; ++i){
		for( j= 0; j< nx; ++j){
			for( ii = i; ii>=0 ; --ii){
				for(jj = j; jj>=0; --jj){
					std::cout<<std::endl;
					sRectc1[j + i*ny] += data[0 + 3 * jj + 3 * nx * ii];
					std::cout<<"value of c1 is "<<data[0 + 3 * jj + 3 * nx * ii]<<"     ";
					sRectc2[j + i*ny] += data[1 + 3 * jj + 3 * nx * ii];
					std::cout<<"value of c2 is "<<data[1 + 3 * jj + 3 * nx * ii]<<"      ";
					sRectc3[j + i*ny] += data[2 + 3 * jj + 3 * nx * ii];
					std::cout<<"value of c3 is "<<data[2 + 3 * jj + 3 * nx * ii]<<"      ";
				}				
			}
			std::cout<<std::endl<<"srectc1="<<sRectc1[j + i*ny]<<"   srectc2="<<sRectc2[j + i*ny]<<"   srectc3="<<sRectc3[j + i*ny];	
		}			
	}
	/*}}}*/	
	std::cout<<std::endl<<"srect1\n";
	for(i = 0; i< ny; ++i){
		for(j = 0; j< nx; ++j){
			std::cout<<sRectc1[j + i*ny]<<" ";
		}
		std::cout<<std::endl;
	}	
	std::cout<<std::endl<<"srect2\n";
	for(i = 0; i< ny; ++i){
                for(j = 0; j< nx; ++j){
                        std::cout<<sRectc2[j + i*ny]<<" ";
                }
                std::cout<<std::endl;
        }
	std::cout<<std::endl<<"srect3\n";
  	for(i = 0; i< ny; ++i){
                for(j = 0; j< nx; ++j){
                        std::cout<<sRectc3[j + i*ny]<<" ";
                }
                std::cout<<std::endl;
        }
	return result;
}
