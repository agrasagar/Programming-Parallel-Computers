#include "cp.h"
#include "cp.h"
#include <iostream>     // std::cout
#include <algorithm>    // std::for_each
#include <vector>       // std::vector
#include <cmath>
#include "../common/vector.h"

void correlate(int ny, int nx, const float* data, float* result) {

	int y, rowStart, rowEnd; 
	int i;
	int rowCount;
	int colCount;
	double meanOfEachRow, meanOfSquareRoot,num;
	std::vector<double> v(nx), vzeroSquaremean(nx);
	int quo = nx / 4;
	int modChck = nx % 4;
	int matSize;
	int travx;	
	if( modChck == 0){
		matSize = quo * ny;
		travx = quo;
	}
	else{
		travx = quo + 1;
		matSize =  ( travx )  * ny ;
	}

	
	double4_t* mat = double4_alloc(matSize);
	rowCount = 0;
	colCount = 0;
//	(void)matSize;

	for(y = 0; y< ny ; ++y){        //traversing through y axis

		rowStart = y*nx;
		rowEnd = rowStart + nx;
		meanOfEachRow = std::accumulate(data+rowStart, data+rowEnd, 0.0)/nx;


		std::transform(data+rowStart, data+rowEnd, v.begin(), [&meanOfEachRow](double val){ return (val - meanOfEachRow); });
		std::transform(v.begin(), v.end(), vzeroSquaremean.begin(), [](double val){ return std::pow(val, 2); });
		meanOfSquareRoot = std::sqrt(std::accumulate(vzeroSquaremean.begin(), vzeroSquaremean.end(), 0.0));
		
		
		if(modChck != 0){
			mat[rowCount + quo] = double4_0;
		}
		for(i = 0; i < nx; ++i) {
			num = v[i]/meanOfSquareRoot;
			if(colCount == 4) {
				colCount = 0;
				rowCount++;
			}
			mat [rowCount][colCount] = num;
			colCount++;	
			
		}
		rowCount++;
		colCount = 0;
	}
	/*	
	//print for mat
	for( i = 0; i < matSize; i++){
		std::cout<<std::endl;
		for(y =0 ; y<4; y++){
			std::cout<<mat[i][y]<<"  ";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
	*/

	//Matrix multiplication
	//matrix nx*ny
	//after multiplication product matrix would be of size ny*ny

	#pragma omp parallel shared (mat, result) 
	{
		int j;
		int k;
		int m;
		#pragma omp for schedule(static, 1)
		for(j = 0; j< ny ; ++j){

			for( k = j; k < ny ; ++k){
				double4_t sum = double4_0;
							
					for(m = 0; m < travx ; ++m){
						
						sum += mat[j*travx + m] * mat[k*travx + m];
					}
				result[k + j*ny] = sum[0] + sum[1] + sum[2] + sum[3];
			}
		}

	}
	free(mat);
}
