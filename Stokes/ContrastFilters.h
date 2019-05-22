#ifndef	CONTRASTE2_H_
#define CONTRASTE2_H_


#include <iostream>
#include <string>
#include <cmath>

#include "Stokes.h"
#include "Polarizers.h"
#include "Window.h"
#include "../Quartic/quartic.h"
#include "LuminanceFilters.h"

using namespace cimg_library;
using namespace std;

unsigned int numNOSOL=0;
float anguloAnterior=0;


//Return the contrast index of the pixel (ix,iy) using the weights and pixels on "pixels"
float contrast(CImg<float> img, vector<Weight> pixels, int ix, int iy){
	float u1=0;
	double total=0;
	for(int i=0;i<pixels.size();i++){
		Weight w=pixels[i];
		int r=w.getRow();
		int c=w.getColumn();
		u1+=w.getWeight()*(luminance(img(r-ix, c-iy, 0, 0), img(r-ix, c-iy, 0, 1), img(r-ix, c-iy, 0, 2)));
		total+=w.getWeight();
	}
	u1=u1*1.0f/pixels.size();
	float k1=0;
	for(int i=0;i<pixels.size();i++){
        Weight w=pixels[i];
        int r=w.getRow();
        int c=w.getColumn();
		k1+=w.getWeight()*pow((luminance(img(r-ix, c-iy, 0, 0), img(r-ix, c-iy, 0, 1),
                               img(r-ix, c-iy, 0, 2)))-u1,2);
	}
	k1=k1*1.0/pixels.size();
	return sqrt(k1);
}

//Return the contrast index of a image
float contrastIndex(CImg<float> img){
	float u1=0;
	for(int i=0;i<img.width();i++){
		for(int j=0;j<img.height();j++){
			u1+=0.2126*img(i,j,0,0)+0.7152*img(i,j,0,1)+0.0722*img(i,j,0,2);
		}
	}
	u1=u1*1.0f/(img.width()*img.height());

	float k1=0;
	for(int i=0;i<img.width();i++){
		for(int j=0;j<img.height();j++){
			k1+=pow((0.2126*img(i,j,0,0)+0.7152*img(i,j,0,1)+0.0722*img(i,j,0,2))-u1,2);
		}
	}
	k1=k1*1.0/(img.width()*img.height());
	return sqrt(k1);
}

/*
    Return the best angle for a linear polarizer to optimize the contrast in a concret pixel.
    (i,j) -> the pixel to optimize,
    maximize -> true to maximize and false to minimize
    win -> Window to process the surrounding pixels
 */
float contrastPixel(StokesImage *vs, int i, int j, bool maximizar, WindowAbs *win){
	//cout << i << " " << j << endl;
	double k1=0;
	double k2=0;
	double k3=0;
	double u1=0,u2=0,u3=0;
	vector<Weight> winWeights = win->getWeights(i, j);
	double total=0;
	int minx=9999,maxx=-9999,miny=9999,maxy=-9999;
	for(int auxV=0;auxV<winWeights.size();auxV++){
		Weight w=winWeights[auxV];
		int r = w.getRow();
		int c = w.getColumn();
		u1+=w.getWeight()*(luminance(vs->atI(r, c, 0), vs->atI(r, c, 1), vs->atI(r, c, 2)));
		u2+=w.getWeight()*(luminance(vs->atQ(r, c, 0), vs->atQ(r, c, 1), vs->atQ(r, c, 2)));
		u3+=w.getWeight()*(luminance(vs->atU(r, c, 0), vs->atU(r, c, 1), vs->atU(r, c, 2)));
		total+=w.getWeight();
		if(r<minx){minx=r;}
		else if(r>maxx){maxx=r;}
		if(c<miny){miny=c;}
		else if(c>maxy){maxy=c;}
	}
	u1=u1/winWeights.size();
	u2=u2/winWeights.size();
	u3=u3/winWeights.size();
	double a=0,b=0,c=0,d=0;
	for(int auxV=0;auxV<winWeights.size();auxV++){
		Weight w=winWeights[auxV];
		int r=w.getRow();
		int col=w.getColumn();
		k1=w.getWeight()*(luminance(vs->atI(r, col, 0), vs->atI(r, col, 1), vs->atI(r, col, 2))-u1);
		k2=w.getWeight()*(luminance(vs->atQ(r, col, 0), vs->atQ(r, col, 1), vs->atQ(r, col, 2))-u2);
		k3=w.getWeight()*(luminance(vs->atU(r, col, 0), vs->atU(r, col, 1), vs->atU(r, col, 2))-u3);
		a+=k1*k3;
		b+=k1*k2;		
		c+=(k3*k3-k2*k2);		
		d+=k2*k3;
	}
	a=a/total;
	b=b/total;
	c=c/total;
	d=d/total;

	float auxa=4*d*d+c*c,auxb=2*a*d+2*c*b,auxc=a*a-4*d*d+b*b-c*c,auxd=-(2*a*d+2*c*b),auxf=d*d-b*b;
	
	std::complex<double> * solutions;
	solutions = solve_quartic(auxb/auxa, auxc/auxa, auxd/auxa, auxf/auxa);
	vector<float> alphas = vector<float>(0);
	for(int k=0;k<4;k++){
		if(solutions[k].imag()==0){
			//REAL ROOT
			alphas.push_back(acos(solutions[k].real())/2);
			alphas.push_back(-acos(solutions[k].real())/2);
		}
	}

	for(int k=0;k<alphas.size();k++){
		if(alphas[k]!=alphas[k]){
			alphas.erase(alphas.begin()+k);
			k--;
		}
	}
	delete[] solutions;
	//NO REAL ROOTS
	if(alphas.size()<=2){
		numNOSOL++;
		vector<float> constAlphas={0,20.0f*cimg::PI/180.0f,40.0f*cimg::PI/180.0f,60.0f*cimg::PI/180.0f,80.0f*cimg::PI/180.0f,100.0f*cimg::PI/180.0f,120.0f*cimg::PI/180.0f,140.0f*cimg::PI/180.0f,160.0f*cimg::PI/180.0f};
		
		for(int k=0;k<constAlphas.size();k++){
			alphas.push_back(constAlphas[k]);
		}
	}

	float max=-99999,min=99999;
	float maxI=0,minI=0;
	float final=0;

	for(int ind=0;ind<alphas.size();ind++){
		CImg<float> filtrada(maxx-minx+1,maxy-miny+1,1,3,0);
		for(int auxi=0;auxi<winWeights.size();auxi++){
			Weight w=winWeights[auxi];
			int r=w.getRow();
			int c=w.getColumn();
			int ai=r-minx, aj=c-miny;
			for(int k=0;k<3;k++){
				filtrada(ai,aj,0,k)=fmin(1,fmax(0,(vs->atI(ai,aj,k)+vs->atQ(ai,aj,k)*cos(2*alphas[ind])+vs->atU(ai,aj,k)*sin(2*alphas[ind]))/2))*255;
			}
		}
		float val= contrast(filtrada, winWeights, minx, miny);
		if(val>max){
			max=val; maxI=ind;
		}
		if(val<min){
			min=val; minI=ind;
		}
	}
	if(maximizar){
		final=alphas[maxI];
	}
	else{
		final=alphas[minI];
	}
	return final;
}

class contrastLocalThread{
	public:
		StokesImage vs;
		mutex mxN;
		CImg<float> resul;
		CImg<float> angles;
		int n, nThreads;
		bool maximize;
		thread * vt[MAX_THREADS];
		WindowAbs * win;
	public:
    contrastLocalThread(StokesImage vs,int nThreads,bool maximize, WindowAbs * win){
			this->vs = vs;
			this->nThreads=nThreads;
			resul= CImg<float>(vs.Stokes.width(),vs.Stokes.height(),1,3,0);
			angles = CImg<float>(vs.Stokes.width(),vs.Stokes.height(),1,3,0);
			n=0;
			this->maximize=maximize;
			this->win=win;
		}

		void filterRow(){
			int i;
			bool maximize;
			mxN.lock();
			i=n;n++;maximize=this->maximize;
			mxN.unlock();
			while(i<vs.Stokes.width()){
				for(int j=0;j<vs.Stokes.height();j++){
					float final = contrastPixel(&vs,i,j,maximize,win);

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

		CImg<float> apply(){
			for(int i=0;i<nThreads;i++){
				vt[i]= new thread( [this] {this->filterRow();});
			}
			for(int i=0;i<nThreads;i++){
				vt[i]->join();
			}

			anglesIMG().save_jpeg("angles.jpg");
			return resul;
		}

		CImg<float> anglesIMG(){
            angles = angles.HSVtoRGB();
			return angles;
		}
};

/*
	Return the image filtered with a different linear polarizer per pixel to maximize or minimize the contrast
    maximize -> true to maximize, false to minimize
    win -> window to use for each pixel
    Generate the image "angles.jpg" where each pixel has a color depending on the angle selected, using the HSV color space
*/
CImg<float> localContrast(StokesImage vs, bool maximize, WindowAbs * win){
    contrastLocalThread fp (vs,20,maximize, win);
	return fp.apply();
}

/*
	Return the image filtered with a linear polarizer to maximize or minimize the contrast
    maximize -> true to maximize, false to minimize
*/
CImg<float> globalContrast(StokesImage vs, bool maximizar){
	double k1=0;
	double k2=0;
	double k3=0;
	double u1=0,u2=0,u3=0;
	vector<Weight> winWeights(0);
	for(int i=0;i<vs.Stokes.width();i++){
		for(int j=0;j<vs.Stokes.height();j++){
			winWeights.push_back(Weight(1,i,j));
		}
	}
	double total=0;
	for(int auxV=0;auxV<winWeights.size();auxV++){
		Weight w=winWeights[auxV];
		//cout << "LOL " << indi << " - " << indj <<endl;
		int r = w.getRow();
		int c = w.getColumn();
		u1+=w.getWeight()*(luminance(vs.atI(r, c, 0), vs.atI(r, c, 1), vs.atI(r, c, 2)));
		u2+=w.getWeight()*(luminance(vs.atQ(r, c, 0), vs.atQ(r, c, 1), vs.atQ(r, c, 2)));
		u3+=w.getWeight()*(luminance(vs.atU(r, c, 0), vs.atU(r, c, 1), vs.atU(r, c, 2)));
		total+=w.getWeight();
	}
	u1=u1/total;
	u2=u2/total;
	u3=u3/total;
	double a=0,b=0,c=0,d=0;
	for(int auxV=0;auxV<winWeights.size();auxV++){
		Weight w=winWeights[auxV];
		int r = w.getRow();
		int col = w.getColumn();
		//cout << "LOL " << indi << " - " << indj <<endl;
		k1=w.getWeight()*((luminance(vs.atI(r, col, 0), vs.atI(r, col, 1), vs.atI(r, col, 2)))-u1);
		k2=w.getWeight()*((luminance(vs.atQ(r, col, 0), vs.atQ(r, col, 1), vs.atQ(r, col, 2)))-u2);
		k3=w.getWeight()*((luminance(vs.atU(r, col, 0), vs.atU(r, col, 1), vs.atU(r, col, 2)))-u3);
		a+=k1*k3;
		b+=k1*k2;		
		c+=(k3*k3-k2*k2);		
		d+=k2*k3;
	}
	a=a/total;
	b=b/total;
	c=c/total;
	d=d/total;
	/*k1=(vs->atI(i,j,0)+vs->atI(i,j,1)+vs->atI(i,j,2))-u1/total;
	k2=(vs->atQ(i,j,0)+vs->atQ(i,j,1)+vs->atQ(i,j,2))-u2/total;
	k3=(vs->atU(i,j,0)+vs->atU(i,j,1)+vs->atU(i,j,2))-u3/total;*/
	float auxa=4*d*d+c*c,auxb=2*a*d+2*c*b,auxc=a*a-4*d*d+b*b-c*c,auxd=-(2*a*d+2*c*b),auxf=d*d-b*b;
	std::complex<double> * solutions;
	solutions = solve_quartic(auxb/auxa, auxc/auxa, auxd/auxa, auxf/auxa);
	vector<float> alphas = vector<float>(0);
	for(int k=0;k<4;k++){
		if(solutions[k].imag()==0){
			//REAL ROOT
			alphas.push_back(acos(solutions[k].real())/2);
			alphas.push_back(-acos(solutions[k].real())/2);
		}
	}

	//Check null values
	for(int k=0;k<alphas.size();k++){
		if(alphas[k]!=alphas[k]){
			alphas.erase(alphas.begin()+k);
			k--;
		}
	}
	delete[] solutions;

	if(alphas.size()<=2){
		cout << "NO REAL ROOTS. SEARCHING IN A SET OF PREDETERMINATED ANGLES" << endl;
		vector<float> constAlphas={0,20.0f*cimg::PI/180.0f,40.0f*cimg::PI/180.0f,60.0f*cimg::PI/180.0f,80.0f*cimg::PI/180.0f,100.0f*cimg::PI/180.0f,120.0f*cimg::PI/180.0f,140.0f*cimg::PI/180.0f,160.0f*cimg::PI/180.0f};
		
		for(int k=0;k<constAlphas.size();k++){
			alphas.push_back(constAlphas[k]);
		}
	}

	float max=-99999,min=99999;
	float maxI=0,minI=0;	
	float final=0;

	for(int ind=0;ind<alphas.size();ind++){
		float val=contrastIndex(linearFilter(vs, alphas[ind] * 180.0f / cimg::PI));
		cout << alphas[ind]*180/cimg::PI << " -> " << val << endl;
		if(val>max){
			max=val; maxI=ind;
		}
		if(val<min){
			min=val; minI=ind;
		}
	}
	if(maximizar){
		final=alphas[maxI];
	}
	else{
		final=alphas[minI];
	}
	cout << "Angle selected = " << final*180/cimg::PI << endl;
	return linearFilter(vs, final * 180 / cimg::PI);
}

#endif