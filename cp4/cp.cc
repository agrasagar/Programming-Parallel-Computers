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
	int quo = nx / 8;
	int modChck = nx % 8;
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

	
	float8_t* mat = float8_alloc(matSize);
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
			mat[rowCount + quo] = float8_0;
		}
		for(i = 0; i < nx; ++i) {
			num = v[i]/meanOfSquareRoot;
			if(colCount == 8) {
				colCount = 0;
				rowCount++;
			}
			mat [rowCount][colCount] = num;
			colCount++;	
			
		}
		rowCount++;
		colCount = 0;
	}

	//Matrix multiplication
	//matrix nx*ny
	//after multiplication product matrix would be of size ny*ny
	
	
	

	#pragma omp parallel  
	{	
		float8_t sjk;
	        float8_t sjkp;
        	float8_t sjkpp;
        	float8_t sjpk;
        	float8_t sjpkp;
        	float8_t sjpkpp;
        	float8_t sjppk;
        	float8_t sjppkp;
        	float8_t sjppkpp;
		int j;
		int k;
		int m;
		bool endChck = false;
		int endVal = ny%3;
		#pragma omp for schedule(static, 1)
		for(j = 0; j< ny ; j=j+3 ){
			if( (j+3)> ny){
				endChck = true;
			}	
			for( k = j; k < ny ; k=k+3 ){
				sjk = float8_0;
                                sjkp = float8_0;
                                sjkpp = float8_0;
                                sjpk = float8_0;
                                sjpkp = float8_0;
                                sjpkpp = float8_0;
                                sjppk = float8_0;
                                sjppkp = float8_0;
                                sjppkpp = float8_0;
				if(endChck != true && (k+3)<ny){
					for(m = 0; m < travx ; ++m){
						sjk += mat[(j)*travx + m] * mat[(k)*travx + m];
						sjkp +=  mat[(j)*travx + m] * mat[(k+1)*travx + m];
						sjkpp += mat[(j)*travx + m] * mat[(k+2)*travx + m];
						
						sjpk +=  mat[(j+1)*travx + m] * mat[(k)*travx + m];
						sjpkp +=  mat[(j+1)*travx + m] * mat[(k+1)*travx + m];
						sjpkpp +=  mat[(j+1)*travx + m] * mat[(k+2)*travx + m];
						
						sjppk +=  mat[(j+2)*travx + m] * mat[(k)*travx + m];
						sjppkp +=  mat[(j+2)*travx + m] * mat[(k+1)*travx + m];
						sjppkpp +=  mat[(j+2)*travx + m] * mat[(k+2)*travx + m];
					}
                        	        result[(k) + (j)*ny] = sjk[0] + sjk[1] + sjk[2] + sjk[3] + sjk[4] + sjk[5] + sjk[6] + sjk[7];
					result[(k+1) + (j)*ny] = sjkp[0] + sjkp[1] + sjkp[2] + sjkp[3] + sjkp[4] + sjkp[5] + sjkp[6] + sjkp[7];
					result[(k+2) + (j)*ny] = sjkpp[0] + sjkpp[1] + sjkpp[2] + sjkpp[3] + sjkpp[4] + sjkpp[5] + sjkpp[6] + sjkpp[7];
				
					result[(k) + (j+1)*ny] = sjpk[0] + sjpk[1] + sjpk[2] + sjpk[3] + sjpk[4] + sjpk[5] + sjpk[6] + sjpk[7];
					result[(k+1) + (j+1)*ny] = sjpkp[0] + sjpkp[1] + sjpkp[2] + sjpkp[3] + sjpkp[4] + sjpkp[5] + sjpkp[6] + sjpkp[7];
					result[(k+2) + (j+1)*ny] = sjpkpp[0] + sjpkpp[1] + sjpkpp[2] + sjpkpp[3] + sjpkpp[4] + sjpkpp[5] + sjpkpp[6] + sjpkpp[7];
					
					result[(k) + (j+2)*ny] = sjppk[0] + sjppk[1] + sjppk[2] + sjppk[3] + sjppk[4] + sjppk[5] + sjppk[6] + sjppk[7];
					result[(k+1) + (j+2)*ny] = sjppkp[0] + sjppkp[1] + sjppkp[2] + sjppkp[3] + sjppkp[4] + sjppkp[5] + sjppkp[6] + sjppkp[7];
					result[(k+2) + (j+2)*ny] = sjppkpp[0] + sjppkpp[1] + sjppkpp[2] + sjppkpp[3] + sjppkpp[4] + sjppkpp[5] + sjppkpp[6] + sjppkpp[7];
				}
				else if(endChck != true && (k+3)>ny){
					if(endVal == 1){
						for(m = 0; m < travx ; ++m){
							sjk += mat[(j)*travx + m] * mat[(k)*travx + m];
							sjpk +=  mat[(j+1)*travx + m] * mat[(k)*travx + m];
							sjppk +=  mat[(j+2)*travx + m] * mat[(k)*travx + m];
						}
						result[(k) + (j)*ny] = sjk[0] + sjk[1] + sjk[2] + sjk[3] + sjk[4] + sjk[5] + sjk[6] + sjk[7];
						result[(k) + (j+1)*ny] = sjpk[0] + sjpk[1] + sjpk[2] + sjpk[3] + sjpk[4] + sjpk[5] + sjpk[6] + sjpk[7];
						result[(k) + (j+2)*ny] = sjppk[0] + sjppk[1] + sjppk[2] + sjppk[3] + sjppk[4] + sjppk[5] + sjppk[6] + sjppk[7];
						break;
					}
					if(endVal == 2){
						for(m = 0; m < travx ; ++m){
                                                        sjk += mat[(j)*travx + m] * mat[(k)*travx + m];
                                                        sjpk +=  mat[(j+1)*travx + m] * mat[(k)*travx + m];
                                                        sjppk +=  mat[(j+2)*travx + m] * mat[(k)*travx + m];
							
							sjkp +=  mat[(j)*travx + m] * mat[(k+1)*travx + m];
							sjpkp +=  mat[(j+1)*travx + m] * mat[(k+1)*travx + m];
							sjppkp +=  mat[(j+2)*travx + m] * mat[(k+1)*travx + m];
							
                                                }
                                                result[(k) + (j)*ny] = sjk[0] + sjk[1] + sjk[2] + sjk[3] + sjk[4] + sjk[5] + sjk[6] + sjk[7];
                                                result[(k) + (j+1)*ny] = sjpk[0] + sjpk[1] + sjpk[2] + sjpk[3] + sjpk[4] + sjpk[5] + sjpk[6] + sjpk[7];
                                                result[(k) + (j+2)*ny] = sjppk[0] + sjppk[1] + sjppk[2] + sjppk[3] + sjppk[4] + sjppk[5] + sjppk[6] + sjppk[7];
						
						result[(k+1) + (j)*ny] = sjkp[0] + sjkp[1] + sjkp[2] + sjkp[3] + sjkp[4] + sjkp[5] + sjkp[6] + sjkp[7];
						result[(k+1) + (j+1)*ny] = sjpkp[0] + sjpkp[1] + sjpkp[2] + sjpkp[3] + sjpkp[4] + sjpkp[5] + sjpkp[6] + sjpkp[7];
						result[(k+1) + (j+2)*ny] = sjppkp[0] + sjppkp[1] + sjppkp[2] + sjppkp[3] + sjppkp[4] + sjppkp[5] + sjppkp[6] + sjppkp[7];
						break;
					}
				}
				else if(endChck == true ){
					if(endVal == 1){
						for(m = 0; m < travx ; ++m){
							sjk += mat[(j)*travx + m] * mat[(k)*travx + m];
						}
						result[(k) + (j)*ny] = sjk[0] + sjk[1] + sjk[2] + sjk[3] + sjk[4] + sjk[5] + sjk[6] + sjk[7];

					}
					if(endVal ==2){
						for(m = 0; m < travx ; ++m){
                                                        sjk += mat[(j)*travx + m] * mat[(k)*travx + m];
                                                        sjkp +=  mat[(j)*travx + m] * mat[(k+1)*travx + m];
							
							sjpk +=  mat[(j+1)*travx + m] * mat[(k)*travx + m];
	                                                sjpkp +=  mat[(j+1)*travx + m] * mat[(k+1)*travx + m];
					        }
                                                result[(k) + (j)*ny] = sjk[0] + sjk[1] + sjk[2] + sjk[3] + sjk[4] + sjk[5] + sjk[6] + sjk[7];
                                                result[(k+1) + (j)*ny] = sjkp[0] + sjkp[1] + sjkp[2] + sjkp[3] + sjkp[4] + sjkp[5] + sjkp[6] + sjkp[7];
						
						result[(k) + (j+1)*ny] = sjpk[0] + sjpk[1] + sjpk[2] + sjpk[3] + sjpk[4] + sjpk[5] + sjpk[6] + sjpk[7];
        	                                result[(k+1) + (j+1)*ny] = sjpkp[0] + sjpkp[1] + sjpkp[2] + sjpkp[3] + sjpkp[4] + sjpkp[5] + sjpkp[6] + sjpkp[7];

					}
				}			

			}//for loop end for k
		}
	}
	free(mat);
}
