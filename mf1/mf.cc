#include "mf.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>

void mf(int ny, int nx, int hy, int hx, const float* in, float* out) {
    int size;
    std::vector<float> medianFilter;
    float median;
    int xy;
   // medianFilter.reserve(hx*hy);   

       // #pragma omp parallel for      
        for(xy = 0; xy < nx*ny; ++xy){   //for traversing through x axis
        int x = xy / ny;
        int y = xy % ny;
	int window_y;
	int window_x;
        medianFilter.clear();

        for(window_y = y - hy; window_y<= y + hy ; window_y++){ //traversing window y axis

          if(window_y >= 0 && window_y < ny){                     //check if window y exceeds given values

            for(window_x = x - hx; window_x<= x + hx; window_x++){  //traversing window x axis

              if(window_x >= 0 && window_x< nx){                  //check if window x exceeds given values

                  medianFilter.push_back(in[window_x + window_y*nx]); //filling up median filter vector

              }
            }

          }

        }

        size = medianFilter.size()/2;
        std::nth_element(medianFilter.begin(), medianFilter.begin() + (size), medianFilter.end());
        median = medianFilter.at(size);
        if ((medianFilter.size() % 2) == 0){
          std::nth_element(medianFilter.begin(), medianFilter.begin()+(size-1), medianFilter.begin()+size);
          median = (median + medianFilter[size-1])/2.0;

        }
        out[x + y*nx] = median;

      }

}

