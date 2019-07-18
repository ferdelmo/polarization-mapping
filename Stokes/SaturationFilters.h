#ifndef SATURACION_H_
#define SATURACION_H_

#include <iostream>
#include <string>
#include <cmath>

#include "Stokes.h"
#include "Polarizers.h"
#include "../Quartic/quartic.h"
#include "LuminanceFilters.h"

using namespace std;


unsigned int numNOSOL2=0; //para contar las veces que Qurtic no encuentra soluciones reales

//Return the saturation value for a color in RGB
float saturation(float R, float G, float B){
	float med=(R+G+B)/3;
	float sat=0;
	float img[]={R,G,B};
	for(int k=0;k<3;k++){
		sat+=(img[k]-med)*(img[k]-med);
	}
	return sqrt(sat/3);
}

float saturation(CImg<float> img, vector<Weight> winWeights, int ix, int iy){
	float resul=0;
	float total=0;
	for(int auxV=0;auxV<winWeights.size();auxV++){
		Weight w = winWeights[auxV];
		int r = w.getRow();
		int c = w.getColumn();
		resul+=w.getWeight()*saturation(img(r-ix,c-iy,0,0),img(r-ix,c-iy,0,1),img(r-ix,c-iy,0,2));
		total+=w.getWeight();
	}
	return resul/total;
}
//return the saturation index of a image (mean of the saturation of its pixels)
float saturationIndex(CImg<float> img){
	float aux=0;
	for(int i=0;i<img.width();i++){
		for(int j=0;j<img.height();j++){
			aux+=saturation(img(i,j,0),img(i,j,1),img(i,j,2));
		}
	}
	return aux/(img.width()*img.height());
}

/*
    Return the best angle for a linear polarizer to optimize the saturation in a concret pixel.
    (i,j) -> the pixel to optimize,
    maximize -> true to maximize and false to minimize
 */
float saturationPixel(StokesImage * vs, int i, int j, bool maximize){
	double iL=0,iQ=0,iU=0;
	for(int k=0;k<3;k++){
		iL+=vs->atI(i,j,k);
		iQ+=vs->atQ(i,j,k);
		iU+=vs->atU(i,j,k);
	}
	iL=iL/3.0f;iQ=iQ/3.0f;iU=iU/3.0f;
	double I=0, Q=0, U=0;

	double a=I*U, b=I*Q, c=(U*U-Q*Q), d=Q*U;
	for(int k=0;k<3;k++){
		I=(vs->atI(i,j,k)-iL);
		Q=(vs->atQ(i,j,k)-iQ);
		U=(vs->atU(i,j,k)-iU);
		a+=(I*U);
		b+=(I*Q);
		c+=(U*U-Q*Q);
		d+=(Q*U);

	}
	double auxa=4*d*d+c*c,auxb=2*a*d+2*c*b,auxc=a*a-4*d*d+b*b-c*c,auxd=-(2*a*d+2*c*b),auxf=d*d-b*b;

	std::complex<double> * solutions;
	solutions = solve_quartic(auxb/auxa, auxc/auxa, auxd/auxa, auxf/auxa);
	vector<float> maximos = vector<float>(0);
	for(int k=0;k<4;k++){
		if(solutions[k].imag()==0){
			//REAL ROOTS
			maximos.push_back(acos(solutions[k].real())/2);
			maximos.push_back(-acos(solutions[k].real())/2);
		}
	}
	//Check nulls values
	for(int k=0;k<maximos.size();k++){
		if(maximos[k]!=maximos[k]){
			//IS NAN
			//cout << "NAAAAN" << endl;
			maximos[k]=0;
		}
	}

	delete[] solutions;
	float max=-9999,min=9999;
	float maxI=0,minI=0;
		
	
	if(maximos.size()<=2){
		//NO REAL SOLUTIONS
		numNOSOL2++;
		int number=9;
		for(int k=0;k<number;k++){
			maximos.push_back(k*180.f/number*cimg::PI/180.0f);
		}
	}

	for(int k=0;k<maximos.size();k++){
		float auxi=2*maximos[k];
		float R=fmin(1,fmax(0,(vs->atI(i,j,0)+vs->atQ(i,j,0)*cos(auxi)+vs->atU(i,j,0)*sin(auxi))/2));
		float G=fmin(1,fmax(0,(vs->atI(i,j,1)+vs->atQ(i,j,1)*cos(auxi)+vs->atU(i,j,1)*sin(auxi))/2));
		float B=fmin(1,fmax(0,(vs->atI(i,j,2)+vs->atQ(i,j,2)*cos(auxi)+vs->atU(i,j,2)*sin(auxi))/2));
		float sat = saturation(R*255,G*255,B*255);
		if(sat > max){
			max=sat;
			maxI=k;
		}
		if(sat < min){
			min=sat;
			minI=k;
		}
	}
	float final=0;
	if(maximize){
		final=maximos[maxI];
	}
	else{
		final=maximos[minI];
	}

	return final;
}

/*
    Return the best angle for a linear polarizer to optimize the saturation in a concret pixel
    and the pixels that surround it, definede by a window.
    (i,j) -> the pixel to optimize,
    maximize -> true to maximize and false to minimize
    win -> window
 */
float saturationWindow(StokesImage * vs, int i, int j, bool maximize, WindowAbs * win){
	vector<Weight> winWeights = win->getWeights(i,j);
	double total=0;
	int minx=9999,maxx=-9999,miny=9999,maxy=-9999;
	for(int auxV=0;auxV<winWeights.size();auxV++){
		Weight w = winWeights[auxV];
		int r = w.getRow();
		int c = w.getColumn();
		/*for(int k=0;k<3;k++){
			iL+=p.p*vs->atI(i,j,k)*255;
			iQ+=p.p*vs->atQ(i,j,k)*255;
			iU+=p.p*vs->atU(i,j,k)*255;
		}*/
		total+=w.getWeight();
		if(r<minx){minx=r;}
		else if(r>maxx){maxx=r;}
		if(c<miny){miny=c;}
		else if(c>maxy){maxy=c;}
	}
	
	float I=0, Q=0, U=0;
	float a=I*U, b=I*Q, c=(U*U-Q*Q), d=Q*U;
	for(int auxV=0;auxV<winWeights.size();auxV++){
		Weight w = winWeights[auxV];
		int r = w.getRow();
		int col = w.getColumn();

		float iL=0,iQ=0,iU=0;
		for(int k=0;k<3;k++){
			iL+=vs->atI(r,col,k)*255;
			iQ+=vs->atQ(r,col,k)*255;
			iU+=vs->atU(r,col,k)*255;
		}
		iL=iL/3.0f;iQ=iQ/3.0f;iU=iU/3.0f;
		for(int k=0;k<3;k++){
			I=(vs->atI(r,col,k)*255-iL);
			Q=(vs->atQ(r,col,k)*255-iQ);
			U=(vs->atU(r,col,k)*255-iU);
			a+=w.getWeight()/total*(I*U);
			b+=w.getWeight()/total*(I*Q);
			c+=w.getWeight()/total*(U*U-Q*Q);
			d+=w.getWeight()/total*(Q*U);

		}
	}
	float auxa=4*d*d+c*c,auxb=2*a*d+2*c*b,auxc=a*a-4*d*d+b*b-c*c,auxd=-(2*a*d+2*c*b),auxf=d*d-b*b;
	
	std::complex<double> * solutions;
	solutions = solve_quartic(auxb/auxa, auxc/auxa, auxd/auxa, auxf/auxa);
	vector<float> maximos = vector<float>(0);
	for(int k=0;k<4;k++){
		if(solutions[k].imag()==0){
			//REAL ROOT
			maximos.push_back(acos(solutions[k].real())/2);
			maximos.push_back(-acos(solutions[k].real())/2);		
		}
	}

	delete[] solutions;
	float max=-9999,min=9999;
	float maxI=0,minI=0;
	for(int k=0;k<maximos.size();k++){
		if(maximos[k]!=maximos[k]){
			maximos[k]=0;
		}
	}
	if(maximos.size()<=2){
		numNOSOL2++;
		vector<float> auxmaximos={0,20.0f*cimg::PI/180.0f,40.0f*cimg::PI/180.0f,60.0f*cimg::PI/180.0f,80.0f*cimg::PI/180.0f,100.0f*cimg::PI/180.0f,120.0f*cimg::PI/180.0f,140.0f*cimg::PI/180.0f,160.0f*cimg::PI/180.0f};
	
		for(int k=0;k<auxmaximos.size();k++){
			maximos.push_back(auxmaximos[k]);
		}
	}
	for(int ind=0;ind<maximos.size();ind++){
		CImg<float> filtrada(maxx-minx+1,maxy-miny+1,1,3,0);
		for(int auxi=0;auxi<winWeights.size();auxi++){
			Weight w=winWeights[auxi];
			int r = w.getRow();
			int c = w.getColumn();
			int ai=r-minx, aj=c-miny;
			for(int k=0;k<3;k++){
				filtrada(ai,aj,0,k)=fmin(1,fmax(0,(vs->atI(r,c,k)+vs->atQ(r,c,k)*cos(2*maximos[ind])+vs->atU(r,c,k)*sin(2*maximos[ind]))/2))*255;
			}
		}

		float sat = saturation(filtrada, winWeights, minx,miny);

		if(sat > max){
			max=sat;
			maxI=ind;
		}
		if(sat < min){
			min=sat;
			minI=ind;
		}
	}
	//cout << "FIN FOR MAXIMOS" << endl;
	float final=0;
	if(maximize){
		final=maximos[maxI];
	}
	else{
		final=maximos[minI];
	}
	return final;
}


/*
	Return the image filtered with a linear polarizer to maximize or minimize the saturation
    maximize -> true to maximize, false to minimize
*/
CImg<float> globalSaturation(StokesImage vs, bool maximize){
	CImg<float> resul(vs.Stokes.width(),vs.Stokes.height(),1,3,200);

	double I=0, Q=0, U=0;
	double a=I*U, b=I*Q, c=(U*U-Q*Q), d=Q*U;
	for(int i=0;i<resul.width();i++){
		//cout << "FILA " << i << endl;
		for(int j=0;j<resul.height();j++){
			double iL=0,iQ=0,iU=0;
			for(int k=0;k<3;k++){
				iL+=vs.atI(i,j,k);
				iQ+=vs.atQ(i,j,k);
				iU+=vs.atU(i,j,k);
			}		
			iL=iL/3.0f;iQ=iQ/3.0f;iU=iU/3.0f;
			for(int k=0;k<3;k++){
				I=255*(vs.atI(i,j,k)-iL);
				Q=255*(vs.atQ(i,j,k)-iQ);
				U=255*(vs.atU(i,j,k)-iU);
				a+=2*(I*U)/(resul.width()*resul.height());
				b+=2*(I*Q)/(resul.width()*resul.height());
				c+=2*(U*U-Q*Q)/(resul.width()*resul.height());
				d+=2*(Q*U)/(resul.width()*resul.height());
			}
		}
	}

	double auxa=4*d*d+c*c,auxb=2*a*d+2*c*b,auxc=a*a-4*d*d+b*b-c*c,auxd=-(2*a*d+2*c*b),auxf=d*d-b*b; 
	std::complex<double> * solutions;
	solutions = solve_quartic(auxb/auxa, auxc/auxa, auxd/auxa, auxf/auxa);
	vector<float> alphas = vector<float>(0);
	for(int k=0;k<4;k++){
		if(solutions[k].imag()==0){
			alphas.push_back(acos(solutions[k].real())/2);
			alphas.push_back(-acos(solutions[k].real())/2);	
		}
	}
	float max=-9999,min=9999;
	float maxI=0,minI=0;
	for(int k=0;k<alphas.size();k++){
		if(alphas[k]!=alphas[k]){
			alphas.erase(alphas.begin()+k);
			k--;
		}
	}
	
	if(alphas.size()<=2){
		cout << "NO REAL ROOTS" << endl;
		int number=9;
		for(int k=0;k<number;k++){
			alphas.push_back(k*180.f/number*cimg::PI/180.0f);
		}
	}

	for(int k=0;k<alphas.size();k++){
		float auxi=alphas[k];
		CImg<float> auxIMG = linearFilter(vs, alphas[k] * 180 / cimg::PI)*255;
		float sat = saturationIndex(auxIMG);
		cout << alphas[k]*180.0f/cimg::PI << " -> " << sat << endl;
		if(sat > max){
			max=sat;
			maxI=k;
		}
		if(sat < min){
			min=sat;
			minI=k;
		}
	}
	float final=0;
	if(alphas.size()>0){
		if(maximize){
			final=alphas[maxI];
		}
		else{
			final=alphas[minI];
		}
	}
	else{
		final=0;
	}
	cout << "Angle selected: " << final*180.0f/cimg::PI << endl;
	for(int i=0;i<resul.width();i++){
		//cout << "FILA " << i << endl;
		for(int j=0;j<resul.height();j++){
			for(int k=0;k<3;k++){
				resul(i,j,k)=fmin(1,fmax(0,(vs.atI(i,j,k)+vs.atQ(i,j,k)*cos(2*final)+vs.atU(i,j,k)*sin(2*final))/2));
			}	
		}
	}
	//nueva.normalize(0,1);
	return resul;
}

class saturationLocalThread{
	public:
		StokesImage vs;
		mutex mxN,mxImg,mxResul;
		CImg<float> resul;
		CImg<float> angles;
		int n, nThreads;
		bool maximize;
		thread * vt[MAX_THREADS];
		WindowAbs * win;
		WindowAbs * winBlur;
		CImg<float> angularMean, normalMean, IMGnormalMean, IMGangularMean;
		bool blur=false;
		bool pixel=true;
	public:
		saturationLocalThread(StokesImage vs,int nThreads,bool maximize, WindowAbs * win, WindowAbs * winBlur){
			this->vs = vs;
			this->nThreads=nThreads;
			resul= CImg<float>(vs.Stokes.width(),vs.Stokes.height(),1,3,0);
			angles = CImg<float>(vs.Stokes.width(),vs.Stokes.height(),1,3,0);
			n=0;
			this->maximize=maximize;
			this->win=win;
			this->pixel=(win==nullptr);
			this->winBlur = winBlur;
			this->blur=!(winBlur==nullptr);

			angularMean = CImg<float>(vs.Stokes.width(),vs.Stokes.height(),1,3,0);
			normalMean = CImg<float>(vs.Stokes.width(),vs.Stokes.height(),1,3,0);
			IMGangularMean = CImg<float>(vs.Stokes.width(),vs.Stokes.height(),1,3,0);
			IMGnormalMean = CImg<float>(vs.Stokes.width(),vs.Stokes.height(),1,3,0);
			//GaussianWindow a(vs.Stokes.width(),vs.Stokes.height(),6,3);
			//winBlur = &a;
			cout << "pixel =" << pixel << " - blur = " << blur << endl; 

		}

		void filterRow(){
			int i;
			bool maximize;
			mxN.lock();
			i=n;n++;maximize=this->maximize;
			mxN.unlock();
			while(i<vs.Stokes.width()){
				for(int j=0;j<vs.Stokes.height();j++){
					float final=0;
					if(pixel){
						final = saturationPixel(&vs,i,j,maximize);
					}
					else{
						final = saturationWindow(&vs,i,j,maximize,win);
					}
					while(final < 0 || final > cimg::PI){
						if(final<0){final+=cimg::PI;}
						else{final-=cimg::PI;}
					}
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

		void blurRow(){
			int i;
			mxN.lock();
			i=n;n++;
			mxN.unlock();
			while(i<vs.Stokes.width()){
				for(int j=0;j<vs.Stokes.height();j++){
					blurAngle(i,j);
					//mxResul.unlock();
				}
				mxN.lock();
				i=n;n++;
				mxN.unlock();
			}
		}
		void blurAngle(int i, int j){
			float angleMed=0;
			float sins = 0;
			float coss = 0;
			float total=0;
			vector<Weight> winWeights = winBlur->getWeights(i, j);
			for(int auxV=0;auxV<winWeights.size();auxV++){
				Weight w = winWeights[auxV];
				int r = w.getRow();
				int c = w.getColumn();
				float rad = angles(r,c,0)/2*cimg::PI/180.0f;
				angleMed+=w.getWeight()*rad; //angles estan entre 0 y 180
				sins += w.getWeight()*sin(2*rad);
				coss += w.getWeight()*cos(2*rad);
				total+=w.getWeight();
			}
			//utilizar una media diferente para angles
			// atan2(sum de los senos, sum de los cosenos)
			angleMed = angleMed/total;
			/*
			normalMean(i,j,0)=2*angleMed*180.0f/cimg::PI;
			normalMean(i,j,1)=1;
			normalMean(i,j,2)=1;
			for(int k=0;k<3;k++){
				IMGnormalMean(i,j,k)=fmin(1,fmax(0,(vs.atI(i,j,k)+vs.atQ(i,j,k)*cos(2*angleMed)+vs.atU(i,j,k)*sin(2*angleMed))/2))*255;
			}
			*/
			float am = atan2(sins/total,coss/total)/2;
			angularMean(i,j,0)=2*am*180.0f/cimg::PI;
			angularMean(i,j,1)=1;
			angularMean(i,j,2)=1;
			for(int k=0;k<3;k++){
				IMGangularMean(i,j,k)=fmin(1,fmax(0,(vs.atI(i,j,k)+vs.atQ(i,j,k)*cos(2*am)+vs.atU(i,j,k)*sin(2*am))/2));
			}
		}
		void blurAngles(){
			cout << "BLUR ANGULOS" << endl;

			resul=resul*255;
			resul.save_jpeg("AngNoMean.jpg");
			n=0;
			for(int i=0;i<nThreads;i++){
				vt[i]= new thread( [this] {this->blurRow();});
			}
			for(int i=0;i<nThreads;i++){
				vt[i]->join();
			}
			angularMean.HSVtoRGB();
			angularMean.save_jpeg("angles.jpg");
			angularMean.RGBtoHSV();
			/*
			normalMean.HSVtoRGB();
			normalMean.save_jpeg("NormalMean.jpg");
			angles.HSVtoRGB();
			IMGnormalMean.save_jpeg("pruebaNormal.jpg");
			IMGangularMean.save_jpeg("pruebaAngular.jpg");
			*/
			//cout << "NO REAL ROOTS: " << numNOSOL2 << "pixels" << endl;
			//cout << "NO REAL ROOTS: " << numNOSOL2*1.0f/(1.0f*angularMean.width()*angularMean.height())*100 << "%" << endl;
			//cout << "SATURATION INDEX: " << saturationIndex(resul) << endl;
			//cout << "INDICE SAT normalMean: " << indiceSaturacion(IMGnormalMean) << endl;
			//cout << "SATURATION INDEX angularMean: " << saturationIndex(IMGangularMean) << endl;

		}

		CImg<float> apply(){
			for(int i=0;i<nThreads;i++){
				vt[i]= new thread( [this] {this->filterRow();});
			}
			for(int i=0;i<nThreads;i++){
				vt[i]->join();
			}

			angles.HSVtoRGB();
			angles.save_jpeg("angles.jpg");
			angles.RGBtoHSV();
			if(blur){
				blurAngles();
				resul = IMGangularMean;
			}
			return resul;
		}
};

/*
	Return the image filtered with a different linear polarizer per pixel to maximize or minimize the saturation
    maximize -> true to maximize, false to minimize
    win -> window to calculate the angle of the pixel, taking account of the pixels that surround it.
    		nullptr to only use the color of the pixel itself
    winBlur -> window use to blur the angles and minimize the noise.
    		nullptr to not blur the angle
    Generate the image "angles.jpg" where each pixel has a color depending on the angle selected, using the HSV color space
*/

CImg<float> localSaturation(StokesImage vs, bool maximize, WindowAbs * win, WindowAbs * winBlur){
	saturationLocalThread fp (vs,20,maximize, win, winBlur);
	CImg<float> resul = fp.apply();
	return resul;
}


#endif