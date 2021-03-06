#include "cp.h"
#include "cp.h"
#include <iostream>     // std::cout
#include <algorithm>    // std::for_each
#include <vector>       // std::vector
#include <cmath>

void correlate(int ny, int nx, const float* data, float* result) {

    int y, rowStart, rowEnd;
    unsigned int i;
    double meanOfEachRow, meanOfSquareRoot,num;
    std::vector<double> v(nx), vzeroSquaremean(nx);
    std::vector <std::vector<double> > matrix(ny);


    for(y = 0; y< ny ; ++y){        //traversing through y axis

      i = 0;
      rowStart = y*nx;
      rowEnd = rowStart + nx;
      meanOfEachRow = std::accumulate(data+rowStart, data+rowEnd, 0.0)/nx;


      std::transform(data+rowStart, data+rowEnd, v.begin(), [&meanOfEachRow](double val){ return (val - meanOfEachRow); });


      std::transform(v.begin(), v.end(), vzeroSquaremean.begin(), [](double val){ return std::pow(val, 2); });


      meanOfSquareRoot = std::sqrt(std::accumulate(vzeroSquaremean.begin(), vzeroSquaremean.end(), 0.0));


    for(i = 0; i < v.size(); ++i) {
      num = v[i]/meanOfSquareRoot;
      matrix[y].push_back(num);
    }
  }
    //Matrix multiplication
    //matrix nx*ny
    //after multiplication product matrix would be of size ny*ny

	#pragma omp parallel shared (matrix, result) 
	{
	int j;
	int k;
	int m;
	#pragma omp for schedule(static, 1)
	for(j = 0; j< ny ; ++j){
				
		for( k = j; k < ny ; ++k){
					double sum;
					sum = 0.0;
					
                        		for(m = 0; m< nx ; ++m){
				
                                		 sum += matrix[j][m] * matrix[k][m];
                        		}	
						result[k + j*ny] = sum;
			}

        }

	}



}
