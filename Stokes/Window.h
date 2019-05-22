#ifndef VENTANA_H_
#define VENTANA_H_


#include <iostream>
#include <vector>
#include <cmath>

#define cimg_display 1
#include "../CImg/CImg.h"

using namespace std;

//class to store a pixel
class Pixel {
public:
    int r,c;
    Pixel(){}
    Pixel(int r, int c){
        this->r = r;
        this->c = c;
    }

    float getRow(){
        return r;
    }

    float getColumn(){
        return c;
    }
};

bool operator==(const Pixel& lhs, const Pixel& rhs){
    return lhs.r == rhs.r && lhs.c == rhs.c;
}

//Store the weight associated to a pixel
class Weight {
	public:
		float w;
        Pixel p;
		Weight(){
			w=1; p=Pixel(0,0);
		}
		Weight(float w, int r, int c){
			this->w=w;
			this->p=Pixel(r,c);
		}

		float getWeight(){
		    return w;
		}

        float getRow(){
            return p.getRow();
        }

        float getColumn(){
            return p.getColumn();
        }
};

//Interface to represent windows
class WindowAbs {
	public:
		int rows, columns, radius;
	public:
        //rows and columns are the size of the target image, used to check the limits of it.
		WindowAbs(int rows, int columns, int radius){
			this->rows=rows;
			this->columns=columns;
			this->radius=radius;
		}
        // Chec if the pixel (r,c) is in the range of the image
		bool inRange(int r, int c){
			return r>=0 && r<rows && c>=0 && c<columns;
		}
		// Return the pixels and their weights of the pixels in the window with center in pixel (r,c)
		virtual vector<Weight> getWeights(int r, int c){
			return {};
		}
};

//Class to represent a Gaussian window
class GaussianWindow : public WindowAbs {
	private:
		int desv;
		vector<vector<float>> pesos;
	public:
        /*
            r and c are the rows and columns of the target image, used to check its limits. Radius is the size
            of the window, and desv the standard deviation of the gaussian
         */
		GaussianWindow(int r,int c, int radius, float desv) : WindowAbs(r,c,radius) {
			this->desv = desv;
			pesos=vector<vector<float>>(radius+1,vector<float>(radius+1,0));
			for(int i=0;i<=radius;i++){
				for(int j=0;j<=radius;j++){
					pesos[i][j]= /*1/(2*cimg::PI*desv*desv)*/exp(-(i*i+j*j)/(2*desv*desv));
					//cout << pesos[i][j] << " ";
				}
				//cout << endl;
			}
		}
        // Return the gaussian weight associated to a pixel x,y (assuming the center is (0,0))
		float gaussianWeight(int x, int y){
			return pesos[abs(x)][abs(y)];
		}
        // Return the pixels and their weights of the pixels in the window with center in pixel (r,c)
		vector<Weight> getWeights(int r, int c){
			vector<Weight> sol = vector<Weight>();
			for(int i=-radius;i<=radius;i++){
				for(int j=-radius;j<=radius;j++){
					if(inRange(r + i, c + j)){
						if(i==0 && j==0){
							sol.push_back(Weight(gaussianWeight(i,j),r,c));
						}
						else{
							sol.push_back(Weight(gaussianWeight(i,j),r+i,c+j));
						}
					}
				}
			}
			return sol;
		}
};

//Represent a basic window where all the weights are 1
class BasicWindow : public WindowAbs {
	public: 
		BasicWindow(int f,int c) : WindowAbs(f,c,1) {

		}

		vector<Weight> getWeights(int r, int c){
			vector<Weight> sol = vector<Weight>();
			for(int i=-1;i<=1;i++){
				for(int j=-1;j<=1;j++){
					if(inRange(r + i, c + j)){
						if(i==0 && j==0){
							sol.push_back(Weight(1,r,c));
						}
						else{
							sol.push_back(Weight(1,r+i,c+j));
						}
					}
				}
			}
			return sol;
		}
};

#endif