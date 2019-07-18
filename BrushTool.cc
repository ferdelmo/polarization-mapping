#include <iostream>

#include "Stokes/Stokes.h"
#include "Stokes/LuminanceFilters.h"
#include "Stokes/ContrastFilters.h"
#include "Stokes/SaturationFilters.h"
#include "Stokes/Polarizers.h"
#include "Stokes/Window.h"


#include <ctime>
#include <string>
#include <vector>
#include <thread>         // std::thread
#include <mutex>          // std::mutex


using namespace cimg_library;
using namespace std;


void color(int x, int y, int R, int G, int B, CImg<float> & img){
	for(int i=x-1;i<x+1;i++){
		for(int j=y-1;j<y+1;j++){
			if(i>=0 && i<img.width() && j>=0 && j<img.height()){
				img(i,j,0)=R;
				img(i,j,1)=G;
				img(i,j,2)=B;
			}
		}
	}
}

void drawcircle(int x0, int y0, int radius, CImg<float> & img, bool maximize, int mode)
{
    int x = radius-1;
    int y = 0;
    int dx = 1;
    int dy = 1;
    int err = dx - (radius << 1);
    int R=255,G=255,B=255;
    if(mode >=0 && mode<=2){
    	if(maximize){
    		G=0;B=0;
    	}
    	else{
    		R=0;G=0;
    	}
    }
    while (x >= y)
    {
        color(x0 + x, y0 + y,R,G,B,img);
        color(x0 + y, y0 + x,R,G,B,img);
        color(x0 - y, y0 + x,R,G,B,img);
        color(x0 - x, y0 + y,R,G,B,img);
        color(x0 - x, y0 - y,R,G,B,img);
        color(x0 - y, y0 - x,R,G,B,img);
        color(x0 + y, y0 - x,R,G,B,img);
        color(x0 + x, y0 - y,R,G,B,img);

        if (err <= 0)
        {
            y++;
            err += dy;
            dy += 2;
        }
        
        if (err > 0)
        {
            x--;
            dx += 2;
            err += dx - (radius << 1);
        }
    }
}

class linealAngulosPara{
	private:
		StokesImage * vs;
		vector<Weight> winWeights;
		CImg<float> angles;
	public:
		void operator()();
		linealAngulosPara(StokesImage * vs, vector<Weight> winWeights, CImg<float> angles){
			this->vs = vs;
			this->winWeights = winWeights;
			this->angles=angles;
		}

		void filtrarLineaN(int i,int n,CImg<float> & l){
			for(int auxV=i;auxV<n;auxV++){
				if(auxV < winWeights.size()){
					Weight w=winWeights[auxV];
					if(w.getWeight()>=0.05f){
						int r=w.getRow();
						int c=w.getColumn();
						float final = angles(r,c,0,0)/2*cimg::PI/180.0f;
						for(int k=0;k<3;k++){
							l(r,c,k)=w.getWeight()*(fmin(1,fmax(0,(vs->atI(r,c,k)+vs->atQ(r,c,k)*cos(2*final)+vs->atU(r,c,k)*sin(2*final))/2))*255)+(1-w.getWeight())*l(r,c,k);
						}
					}
				}
			}
		}

		thread tothread(int i, int n, CImg<float> & l){
			return thread(&linealAngulosPara::filtrarLineaN,this,i,n,std::ref(l));
		}
};

void linealAngulos(StokesImage * vs, int x, int y, WindowAbs * ven, CImg<float> & anterior, CImg<float> & angulos){
	vector<Weight> winWeights = ven->getWeights(x, y);
	linealAngulosPara fp (vs, winWeights,angulos);
	int total = 4;
	thread vt[total];
	//cout << "THREADS: " << total << " - PIXELES: " << winWeights.size()/total << endl;
	for(int i=0;i<total;i++){
		vt[i]= fp.tothread(i*winWeights.size()/total, (i+1)*winWeights.size()/total, anterior);
	}
	for(int i=0;i<total;i++){
		vt[i].join();
	}
}

int main(int argc, char* argv[]) {
	int height = 1080;
	string dir="";
	for(int i=0;i<argc;i++){
		string arg = argv[i];
		if(arg=="-i"){
			i++;
			dir=argv[i];
		}
	}
	StokesImage vs(dir);
	//StokesImage vs(CImg<float>("prueba8/salida/I.jpg"),CImg<float>("prueba8/salida/Q.jpg"),CImg<float>("prueba8/salida/U.jpg"),CImg<float>("prueba8/salida/V.jpg"));
	
	vs.resize(vs.Stokes.width()*1.0f/vs.Stokes.height()*height,height);
	CImg<float> nueva(vs.Stokes.width(),vs.Stokes.height(),1,3,0);

	CImg<float> angBmax, angBmin, angCmax, angCmin, angSmax, angSmin;
	//ANGULOS BRILLO
	/*cout << "CALCULATING LUMINANCE ANGLES" << endl;
	luminanceLocalThreads fpbmax (vs,10,true);
    fpbmax.apply();
	angBmax = fpbmax.angles;

	luminanceLocalThreads fpbmin (vs,10,false);
    fpbmin.apply();
	angBmin = fpbmin.angles;*/
	//ANGULOS CONTRASTE
	cout << "CALCULATING CONTRAST ANGLES" << endl;
	GaussianWindow ac(vs.Stokes.width(),vs.Stokes.height(),10,5);
	WindowAbs * ptrc = &ac;
	contrastLocalThread fpcmax (vs,10,true,ptrc);
	fpcmax.apply();
	angCmax = fpcmax.angles;

	contrastLocalThread fpcmin (vs,10,false,ptrc);
	fpcmin.apply();
	angCmin = fpcmin.angles;
	//ANGULOS SATURACION
	cout << "CALCULATING SATURATION ANGLES" << endl;
	GaussianWindow a2(vs.Stokes.width(),vs.Stokes.height(),15,7.5);
	WindowAbs * ptr2 = &a2;
	saturationLocalThread fpsmax (vs,10,true,nullptr,ptr2);
	fpsmax.apply();
	angSmax = fpsmax.angularMean;
	saturationLocalThread fpsmin (vs,10,false,nullptr,ptr2);
	fpsmin.apply();
	angSmin = fpsmin.angularMean;


	for(int i=0;i<nueva.width();i++){
		for(int j=0;j<nueva.height();j++){
			for(int k=0;k<3;k++){
				nueva(i,j,k)=vs.atI(i,j,k)/2*255;
			}
		}
	}
	CImgDisplay main_disp(nueva,"BRUSH TOOL",0);
	int tam=20;
	bool maximizar=false;
	bool pulsadaM=true;
	int modo=0; //0 brillo, 1 contraste, 2 sat, 3 original, 4 lineal
	bool cambioCon=false;
	bool pulsada = true;
	float angulo=0;
	GaussianWindow a(nueva.width(),nueva.height(),tam,tam/3.5f);
	while (!main_disp.is_closed()){
		clock_t begin = clock();
		CImg<float> mostrar(nueva);
		main_disp.wait();
		clock_t end = clock();
		double msi = double(end-begin)/(CLOCKS_PER_SEC/1000);
		//cout << (main_disp.button()&1) << endl;
		begin = clock();
		if(main_disp.is_keyPADADD()){
			tam++;
			cout << "WINDOW ENLARGED" << endl;
			cout << tam << endl;
			a = GaussianWindow(nueva.width(),nueva.height(),tam,tam/2.5f);
		}
		if(main_disp.is_keyPADSUB()){
			if(tam>2){
				tam--;
				cout << "SHRUNKEN WINDOW" << endl;
				a = GaussianWindow(nueva.width(),nueva.height(),tam,tam/2.5f);
			}
		}
		//MAXIMIZAR O MINIMIZAR
		if(main_disp.is_keyM()){
			if(pulsadaM){
				cambioCon=true;
				maximizar = !maximizar;
				if(maximizar){
					cout << "MAXIMIZING" << endl;
				}
				else{
					cout << "MINIMIZING" << endl;
				}
				pulsadaM=false;
			}
		}
		else{
			pulsadaM=true;
		}

		if(main_disp.is_keyL()){
			if(pulsada){
				cout << "OPTIMIZING LUMINANCE" << endl;
				modo=0;
				pulsada = false;
			}
		}
		else if(main_disp.is_keyC()){
			if(pulsada){
				cout << "OPTIMIZING CONTRAST" << endl;
				modo=1;
				pulsada = false;
			}
		}
		else if(main_disp.is_keyS()){
			if(pulsada){
				cout << "OPTIMIZING SATURATION" << endl;
				modo=2;
				pulsada = false;
			}
		}
		else if(main_disp.is_keyO()){
			if(pulsada){
				cout << "APLICANANDO LA IMAGEN ORIGINAL" << endl;
				modo=3;
				pulsada = false;
			}
		}
		else if(main_disp.is_keyP()){
			if(pulsada){
				cout << "LINEAR POLARIZER. Introduce angle in degrees" << endl;
				string ang="";
				cin >> ang;
				angulo=stof(ang);
				modo=4;
				pulsada = false;
			}
		}
		else if(main_disp.is_keyG()){
			if(pulsada){
				cout << "Path to save the image: " << endl;
				string ruta="";
				cin >> ruta;
				mostrar.save_jpeg(ruta.c_str());
				cout << "IMAGEN GUARDADA" << endl;
				pulsada = false;
			}
		}
		else{
			pulsada = true;
		}
		end = clock();
		double msb = double(end-begin)/(CLOCKS_PER_SEC/1000);
		begin = clock();
		drawcircle(main_disp.mouse_x(),main_disp.mouse_y(),tam,mostrar,maximizar,modo);
		if(modo==0) {
			if(main_disp.button()&1 && main_disp.mouse_y()>=0 && main_disp.mouse_x()>=0){
				int x = main_disp.mouse_x(), y = main_disp.mouse_y();
				luminanceBrush(&vs, maximizar, a.radius,a.getWeights(x,y), nueva);
				/*if(maximizar){
					linealAngulos(&vs,x,y,&a,nueva,angBmax);
				}
				else{
					linealAngulos(&vs,x,y,&a,nueva,angBmin);
				}*/
			}	
		}
		else if(modo==1){
			if(main_disp.button()&1 && main_disp.mouse_y()>=0 && main_disp.mouse_x()>=0){
				int x = main_disp.mouse_x(), y = main_disp.mouse_y();
				if(maximizar){
					linealAngulos(&vs,x,y,&a,nueva,angCmax);
				}
				else{
					linealAngulos(&vs,x,y,&a,nueva,angCmin);
				}
			}
		}
		else if(modo==2){
			if(main_disp.button()&1 && main_disp.mouse_y()>=0 && main_disp.mouse_x()>=0){
				int x = main_disp.mouse_x(), y = main_disp.mouse_y();
				if(maximizar){
					linealAngulos(&vs,x,y,&a,nueva,angSmax);
				}
				else{
					linealAngulos(&vs,x,y,&a,nueva,angSmin);
				}
			}
		}
		else if(modo==3){
			if(main_disp.button()&1 && main_disp.mouse_y()>=0 && main_disp.mouse_x()>=0){
				int x = main_disp.mouse_x(), y = main_disp.mouse_y();
				original(vs,x,y,&a,nueva);
			}
		}
		else if(modo==4){
			if(main_disp.button()&1 && main_disp.mouse_y()>=0 && main_disp.mouse_x()>=0){
				int x = main_disp.mouse_x(), y = main_disp.mouse_y();
                linearWindow(vs, x, y, &a, nueva, angulo);
			}
		}
		main_disp.display(mostrar);
	}
}