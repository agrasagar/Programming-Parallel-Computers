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
    int y, rowStart, rowEnd, i,x,k;
    double meanOfEachRow, meanOfSquareRoot;
    std::vector<double> v(nx), vzeroSquaremean(nx), vrowfinal(nx);
    std::vector <std::vector<double> > matrix(ny);
    std::vector <std::vector<double> > matrixT(nx);


    for(y = 0; y< ny ; ++y){        //traversing through y axis

      i = 0;
      rowStart = y*nx;
      rowEnd = rowStart + nx;
      meanOfEachRow = std::accumulate(data+rowStart, data+rowEnd, 0.0)/nx;

      //std::cout<<std::endl<<meanOfEachRow<<std::endl;

      std::transform(data+rowStart, data+rowEnd, v.begin(), [&meanOfEachRow](double val){ return (val - meanOfEachRow); });
      //std::cout<<v.size()<<std::endl;

      std::transform(v.begin(), v.end(), vzeroSquaremean.begin(), [](double val){ return std::pow(val, 2); });
      //std::cout<<vzeroSquaremean.size()<<std::endl;

      meanOfSquareRoot = std::sqrt(std::accumulate(vzeroSquaremean.begin(), vzeroSquaremean.end(), 0.0));
      //std::cout<<std::endl<<meanOfSquareRoot<<std::endl;

      std::transform(v.begin(), v.end(), vrowfinal.begin(), [&meanOfSquareRoot,&i,&matrixT,&matrix,&y](double val){
        double num;
        num = val/meanOfSquareRoot;
        matrixT[i].push_back(num);
        i++;
        matrix[y].push_back(num);
        return num; });


      //std::cout<<vrowfinal.size()<<std::endl;

    }

    //Matrix multiplication
    //matrix nx*ny
    //matrixT ny*nx
    //after multiplication product matrix would be of size ny*ny

  

       for(y = 0; y< ny ; ++y){
         for(k = 0; k<nx; ++k){
	   for(x = 0; x< ny ; ++x){
		result[x + y*ny] += matrix[y][k]*matrixT[k][x];
           }

         }
       }
}

