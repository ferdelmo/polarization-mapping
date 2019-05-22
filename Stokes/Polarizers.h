#ifndef	LINEARPOLARIZATION_H_
#define LINEARPOLARIZATION_H_


#include <iostream>
#include <string>
#include <cmath>

#include "Stokes.h"
#include "Window.h"

using namespace cimg_library;
using namespace std;

// Change from degrees to radians
float deg2rad(float a){
	return a*cimg::PI/180.0f;
}

/*
	Return an image filtered with a linear polarizer at g degrees
*/
CImg<float> linearFilter(StokesImage vs, float g){
	CImg<float> resul(vs.Stokes.width(),vs.Stokes.height(),1,3,0);
	float rad=deg2rad(g);
	for(int i=0;i<resul.width();i++){
		for(int j=0;j<resul.height();j++){
			for(int k=0;k<3;k++){
                resul(i,j,k)=fmin(1,fmax(0,(vs.atI(i,j,k)+vs.atQ(i,j,k)*cos(2*rad)+vs.atU(i,j,k)*sin(2*rad))/2));
			}
		}
	}
	return resul;
}

/*
	Return an image filtered with an RCP filter
*/
CImg<float> RCP(StokesImage vs){
	CImg<float> resul(vs.Stokes.width(),vs.Stokes.height(),1,3,0);
	for(int i=0;i<resul.width();i++){
		for(int j=0;j<resul.height();j++){
			for(int k=0;k<3;k++){
                resul(i,j,k)=fmin(1,fmax(0,(vs.atI(i,j,k)-vs.atV(i,j,k))/2));
			}
		}
	}
	return resul;
}

/*
	Return an image filtered with an LCP filter
*/
CImg<float> LCP(StokesImage vs){
	CImg<float> nueva(vs.Stokes.width(),vs.Stokes.height(),1,3,0);
	for(int i=0;i<nueva.width();i++){
		for(int j=0;j<nueva.height();j++){
			for(int k=0;k<3;k++){
				nueva(i,j,k)=fmin(1,fmax(0,(vs.atI(i,j,k)+vs.atV(i,j,k))/2));
			}
		}
	}
	return nueva;
}
/*
	Return an image filtered with a custom filter with 3 coeficients for the parameters Q,U and V (c1,c2,c3)
*/
CImg<float> customFilter(StokesImage vs, float c1, float c2, float c3){
	CImg<float> resul(vs.Stokes.width(),vs.Stokes.height(),1,3,0);
	for(int i=0;i<resul.width();i++){
		for(int j=0;j<resul.height();j++){
			for(int k=0;k<3;k++){
                resul(i,j,k)=fmin(1,fmax(0,(vs.atI(i,j,k)+vs.atQ(i,j,k)*c1+vs.atU(i,j,k)*c2+vs.atV(i,j,k)*c3)/2));
			}
		}
	}
	return resul;
}

/*
    Restore the image without any polarizer in the pixel (x,y) and the window win, and store the result in "last"
 */

void original(StokesImage vs, int x, int y, WindowAbs * win, CImg<float> & last){

	vector<Weight> winWeights = win->getWeights(x, y);
	for(int auxV=0;auxV<winWeights.size();auxV++){

		Weight w=winWeights[auxV];
		for(int k=0;k<3;k++){
		    int r = w.getRow();
		    int c = w.getColumn();
			if(w.getWeight()>=0.05f){
                last(r,c,k)=w.getWeight()*(fmin(1,fmax(0,(vs.atI(r,c,k))/2))*255)+(1-w.getWeight())*last(r,c,k);
			}
		}
	}
	for(int k=0;k<3;k++){
        last(x,y,k)=fmin(1,fmax(0,(vs.atI(x,y,k))/2))*255;
	}
}

/*
	Apply a linear polarizer in the pixels of the window "win" with center in pixel (x,y). The changes are stored in
    the image "last". "ang" is the angle of the linear polarizer
*/
void linearWindow(StokesImage vs, int x, int y, WindowAbs *win, CImg<float> &last, float ang){
	vector<Weight> winWeights = win->getWeights(x, y);
    float final = ang;
	for(int auxV=0;auxV<winWeights.size();auxV++){
        Weight w=winWeights[auxV];
		for(int k=0;k<3;k++){
            int r = w.getRow();
            int c = w.getColumn();
			if(w.getWeight()>=0.05f){
				last(r,c,k)=w.getWeight()*(fmin(1,fmax(0,(vs.atI(r,c,k)+vs.atQ(r,c,k)*cos(2*final)+vs.atU(r,c,k)*sin(2*final))/2))*255)+(1-w.getWeight())*last(r,c,k);
			}
		}
	}
	for(int k=0;k<3;k++){
		last(x,y,k)=fmin(1,fmax(0,(vs.atI(x,y,k)+vs.atQ(x,y,k)*cos(2*final)+vs.atU(x,y,k)*sin(2*final))/2))*255;
	}
}


/*
	Return an image filtered with a polarizer composed by a linear polarizer at alpha radians and a shifter at delta radians
*/

CImg<float> circularFilter(StokesImage vs, float alpha, float delta){
	CImg<float> resul(vs.Stokes.width(),vs.Stokes.height(),1,3,0);
	float x=2*alpha;
	float y=delta;
	for(int i=0;i<resul.width();i++){
		for(int j=0;j<resul.height();j++){
			float cQ=(cos(x)-sin(x)*cos(y))/sqrt(2);
			float cU=(cos(x)*cos(y)+sin(x))/sqrt(2);
			float cV=-sin(y)/sqrt(2);
			for(int k=0;k<3;k++){
                resul(i,j,k)=fmin(1,fmax(0,(vs.atI(i,j,k)+cQ*vs.atQ(i,j,k)+cU*vs.atU(i,j,k)+cV*vs.atV(i,j,k))/2));
			}
		}
	}
	return resul;
}

#endif