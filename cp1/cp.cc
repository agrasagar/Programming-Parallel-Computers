#include "cp.h"
#include "cp.h"
#include <iostream>     // std::cout
#include <algorithm>    // std::for_each
#include <vector>       // std::vector
#include <cmath>

void correlate(int ny, int nx, const float* data, float* result) {
    /*
    // FIXME
    for (int i = 0; i < ny * ny; ++i) {
        result[i] = 0.0f;
    }
    */
    int y, rowStart, rowEnd, i,x;
    double meanOfEachRow, meanOfSquareRoot;
    std::vector<double> v(nx), vzeroSquaremean(nx), vrowfinal(nx);
    std::vector <std::vector<double> > matrix(ny, std::vector<double>(nx));
    //std::vector <std::vector<double> > matrixT(nx, std::vector<double>(ny));
    std::vector <std::vector<double> > matrixT(nx);
//    v.reserve(nx);
//    vzeroSquaremean.reserve(nx);
//    vrowfinal.reserve(nx); 
//    matrix.reserve(ny);
//    matrixT.reserve(nx);
    //std::cout<<std::endl<<v.size()<<std::endl;
    for(y = 0; y< ny ; ++y){        //traversing through y axis
      //vzeroSquaremean.clear();
      //vrowfinal.clear();
      i = 0;
      rowStart = y*nx;
      rowEnd = rowStart + nx;
      meanOfEachRow = std::accumulate(data+rowStart, data+rowEnd, 0.0)/nx;
      //std::cout<<std::endl<<meanOfEachRow<<std::endl;

      std::transform(data+rowStart, data+rowEnd, v.begin(), [&meanOfEachRow](double val){ return (val - meanOfEachRow); });
      //std::cout<<v.size()<<std::endl;

      std::transform(v.begin(), v.end(), vzeroSquaremean.begin(), [](double val){ return pow(val, 2); });
      //std::cout<<vzeroSquaremean.size()<<std::endl;

      meanOfSquareRoot = sqrt(std::accumulate(vzeroSquaremean.begin(), vzeroSquaremean.end(), 0.0));
      //std::cout<<std::endl<<meanOfSquareRoot<<std::endl;

      std::transform(v.begin(), v.end(), vrowfinal.begin(), [&meanOfSquareRoot](double val){return val/meanOfSquareRoot; });

      std::copy(vrowfinal.begin(), vrowfinal.end(),matrix[y].begin());

      //std::cout<<vrowfinal.size()<<std::endl;

      //Transpose
      i = 0;
      for(auto& num: vrowfinal){
        matrixT[i].push_back(num);
        i++;
      }

    }

    //Matrix multiplication
    //matrix nx*ny
    //matrixT ny*nx
    //after multiplication product matrix would be of size ny*ny
//	std::cout<<"reach till here"<<std::endl;
  //      std::cout<<matrix.size()<<std::endl;
    //    std::cout<<matrix[0].size()<<std::endl;

//	std::cout<<matrixT.size()<<std::endl;
//	std::cout<<matrixT[0].size()<<std::endl;


       int k;
       for(y = 0; y< ny ; ++y){
         for(k = 0; k<nx; ++k){
	   for(x = 0; x< ny ; ++x){
		result[x + y*ny] += matrix[y][k]*matrixT[k][x];
           }

         }
       }
}

