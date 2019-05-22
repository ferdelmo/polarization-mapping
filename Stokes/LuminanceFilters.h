#ifndef	BRILLO_H_
#define BRILLO_H_


#include <iostream>
#include <string>
#include <cmath>
#include <thread>         // std::thread
#include <mutex>          // std::mutex

#include <valarray>

#include "Stokes.h"
#include "Polarizers.h"
#include "Window.h"

using namespace cimg_library;
using namespace std;

// Return the luminance of a color
float luminance(float R, float G, float B){
	return R*0.2126f+G*0.7152f+B*0.0722f;
}

//Return the luminance index of an image
float luminanceIndex(CImg<float> img){
	float resul=0;
	for(int i=0;i<img.width();i++){
		for(int j=0;j<img.height();j++){
			resul+= luminance(img(i, j, 0, 0), img(i, j, 0, 1), img(i, j, 0, 2));
		}
	}
	return resul/(img.width()*img.height());
}

/*
    Return the best angle for a linear polarizer to optimize the luminance in a concret pixel.
    (i,j) -> the pixel to optimize,
    maximize -> true to maximize and false to minimize
    luminance -> the function that describes luminance
 */

float luminancePixel(StokesImage *vs, int i, int j, bool maximize, float (*luminance)(float R, float G, float B)){
	float k1=luminance(vs->atQ(i,j,0),vs->atQ(i,j,1),vs->atQ(i,j,2));
	float k2=luminance(vs->atU(i,j,0),vs->atU(i,j,1),vs->atU(i,j,2));
	float x1=0;
	float x2=0;
	if(k1!=0){
		x1=atan(k2/k1)/2;
		x2=(atan(k2/k1)+cimg::PI)/2;
	}
	else{
		x1=cimg::PI/2/2;
		x2=(cimg::PI/2+cimg::PI)/2;
	}
	float final=0;
	if(isnan(x1) || isnan(x2)){
		x1=cimg::PI/2/2;
		x2=(cimg::PI/2+cimg::PI)/2;
	}
	if(-k1*cos(2*x1)-k2*sin(2*x1)<0){
		if(maximize){
			final=x1;
		}
		else{
			final=x2;
		}
	}
	else{
		if(maximize){
			final=x2;
		}
		else{
			final=x1;
		}
	}
	return final;
}


/*
	Return the image filtered with a linear polarizer to maximize or minimize the luminance
    maximize -> true to maximize, false to minimize
*/
CImg<float> globalLuminance(StokesImage vs, bool maximize){
	float k1=0;
	float k2=0;
	int tam=vs.Stokes.width()*vs.Stokes.height();
	for(int i=0;i<vs.Stokes.width();i++){
		for(int j=0;j<vs.Stokes.height();j++){
			k1+=(vs.atQ(i,j,0)+vs.atQ(i,j,1)+vs.atQ(i,j,2))/tam;
			k2+=(vs.atU(i,j,0)+vs.atU(i,j,1)+vs.atU(i,j,2))/tam;
		}
	}
	float x1=0;
	float x2=0;
	if(k1!=0){
		x1=atan(k2/k1)/2;
		x2=(atan(k2/k1)+cimg::PI)/2;
	}
	else{
		x1=cimg::PI/2/2;
		x2=(cimg::PI/2+cimg::PI)/2;
	}
	float final=0;

	if(isnan(x1) || isnan(x2)){
		x1=cimg::PI/2/2;
		x2=(cimg::PI/2+cimg::PI)/2;
	}
	if(-k1*cos(2*x1)-k2*sin(2*x1)<0){
		if(maximize){
			final=x1;
		}
		else{
			final=x2;
		}
	}
	else{
		if(maximize){
			final=x2;
		}
		else{
			final=x1;
		}
	}
	cout << "Angle selected = " << final*180/cimg::PI << "º" << endl;
	return linearFilter(vs, final * 180 / cimg::PI);
}

const int MAX_THREADS=200;

//Class to apply the local luminance filter with paralelization
class luminanceLocalThreads{
	public:
		StokesImage vs;
		mutex mxN;
		CImg<float> resul;
		int n, nThreads;
		bool maximize;
		thread * vt[MAX_THREADS];
		CImg<float> angles;
	public:
		luminanceLocalThreads(StokesImage vs,int nThreads,bool maximize){
			this->vs = vs;
			this->nThreads=nThreads;
			resul= CImg<float>(vs.Stokes.width(),vs.Stokes.height(),1,3,0);
			angles= CImg<float>(vs.Stokes.width(),vs.Stokes.height(),1,3,0);
			n=0;
			this->maximize=maximize;
		}

		//process the row n
		void filterRow(){
			int i;
			bool maximizar;
			mxN.lock();
			i=n;n++;maximizar=this->maximize;
			mxN.unlock();
			while(i<vs.Stokes.width()){
				for(int j=0;j<vs.Stokes.height();j++){
					float final = luminancePixel(&vs, i, j, maximizar, &luminance);
					angles(i,j,0)=2*final*180.0f/cimg::PI;
					angles(i,j,1)=1;
					angles(i,j,2)=1;
					//mxResul.lock();
					for(int k=0;k<3;k++){
						resul(i,j,k)=fmin(1,fmax(0,(vs.atI(i,j,k)+vs.atQ(i,j,k)*cos(2*final)+vs.atU(i,j,k)*sin(2*final))/2));
					}
					//mxResul.unlock();
				}
				mxN.lock();
				i=n;n++;
				mxN.unlock();
			}
		}

		// Apply the luminance local filter using nThreads
		CImg<float> apply(){
			for(int i=0;i<nThreads;i++){
				vt[i]= new thread( [this] { this->filterRow();});
			}
			for(int i=0;i<nThreads;i++){
				vt[i]->join();
			}
            anglesIMG().save_jpeg("angles.jpg");
			angles.RGBtoHSV();
			return resul;
		}

		CImg<float> anglesIMG(){
			return angles.HSVtoRGB();
		}
};

/*
	Return the image filtered with a different linear polarizer per pixel to maximize or minimize the luminance
    maximize -> true to maximize, false to minimize
    Generate the image "angles.jpg" where each pixel has a color depending on the angle selected, using the HSV color space
*/
CImg<float> localLuminance(StokesImage vs, bool maximize){
	luminanceLocalThreads fp (vs,25,maximize);
	return fp.apply();
}

//Class to apply the local luminance filter in a window with paralelization (for the brush-tool)
class luminanceFilterWindow{
	private:
		StokesImage * vs;
		bool maximize;
		vector<Weight> winWeights;
	public:
		void operator()();

		luminanceFilterWindow(StokesImage * vs,bool maximize, vector<Weight> winWeights){
			this->vs = vs;
			this->maximize=maximize;
			this->winWeights = winWeights;
		}

		void filterRows(int i, int n, CImg<float> &last){
			for(int auxV=i;auxV<n;auxV++){
				if(auxV < winWeights.size()){
                    Weight w=winWeights[auxV];
                    int r = w.getRow();
                    int c = w.getColumn();

					if(w.getWeight()>=0.05f){
						float final = luminancePixel(vs, r, c, maximize, &luminance);
						for(int k=0;k<3;k++){
							last(r,c,k)=w.getWeight()*(fmin(1,fmax(0,(vs->atI(r,c,k)+vs->atQ(r,c,k)*cos(2*final)+
							vs->atU(r,c,k)*sin(2*final))/2))*255)+(1-w.getWeight())*last(r,c,k);
						}
					}
				}
			}
		}

		thread tothread(int i, int n, CImg<float> & anterior){
			return thread(&luminanceFilterWindow::filterRows,this,i,n,std::ref(anterior));
		}
};
void filtrosBrilloParaleloPincel(StokesImage * vs, bool maximize, int radius, vector<Weight> winWeights, CImg<float> & last){
	luminanceFilterWindow fp (vs,maximize, winWeights);
	int total = radius/2;
	thread vt[total];
	//cout << "THREADS: " << total << " - PIXELES: " << winWeights.size()/total << endl;
	for(int i=0;i<total;i++){
		vt[i]= fp.tothread(i*winWeights.size()/total, (i+1)*winWeights.size()/total, last);
	}
	for(int i=0;i<total;i++){
		vt[i].join();
	}
}


#endif