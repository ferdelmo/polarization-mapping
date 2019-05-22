#include <iostream>

#include "Stokes.h"
#include "LuminanceFilters.h"
#include "ContrastFilters.h"
#include "SaturationFilters.h"
#include "Polarizers.h"
#include "Window.h"


#include <ctime>
#include <string>
#include <vector>
#include <thread>         // std::thread
#include <mutex>          // std::mutex


using namespace cimg_library;
using namespace std;


/*
	Pinta 3 pixeles
*/
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

class Rectangulo{
	private:
		int R=255,G=255,B=255;
	public:
		Pixel i,p;
		bool actu=false;
		Rectangulo(){}
		void dibujar(CImg<float> & img){
			int izq = i.c, der = p.c;
			if(izq>der){
				izq = p.c;
				der = i.c;
			}

			int top = i.f, bot = p.f;
			if(bot>top){
				top = p.f;
				bot = i.f;
			}

			//HORIZONTALES
			for(int aux=izq;aux<der;aux++){
				color(aux,bot,R,G,B,img);
				color(aux,top,R,G,B,img);
			}

			//VERTICALES
			for(int aux=bot;aux<top;aux++){
				color(izq,aux,R,G,B,img);
				color(der,aux,R,G,B,img);
			}
		}

		vector<Pixel> getPixeles(){
			vector<Pixel> resul = {};
			int izq = i.c, der = p.c;
			if(izq>der){
				izq = p.c;
				der = i.c;
			}
			int top = i.f, bot = p.f;
			if(bot>top){
				top = p.f;
				bot = i.f;
			}
			for(int i=izq;i<der;i++){
				for(int j=bot;j<top;j++){
					resul.push_back(Pixel(i,j));
				}
			}
			return resul;
		}
};

/*
	Clase para contar los FPS
*/
class FPS{
	int frames;
	clock_t begin,end;
	bool emp=true,stop = false;
	mutex m;
	thread c;
	double fps;


public:
	FPS(){
		frames=0;
	}
	void empezar(){
		c = thread( [this] {this->contador();});

	}
	void parar(){
		m.lock();
		stop = true;
		m.unlock();
		c.join();
	}
	void contador(){
		m.lock();
		while(!stop){
			m.unlock();
			if(emp){
				begin = clock();
				m.lock();
				frames=0;
				m.unlock();
				emp=false;
			}
			else{
				end = clock();
				double ms = double(end - begin) / (CLOCKS_PER_SEC/1000);
				if(ms>1000){
					m.lock();
					fps = frames * 1.0 / ms *1000;
					m.unlock();
					emp=true;
				}
			}
		}
		m.unlock();
	}
	void showFPS(){
		m.lock();
		cout << "FPS: " << fps << endl;
		m.unlock();
	}
	void frame(){
		m.lock();
		frames++;
		m.unlock();
	}
};

void drawcircle(int x0, int y0, int radius, CImg<float> & img)
{
    int x = radius-1;
    int y = 0;
    int dx = 1;
    int dy = 1;
    int err = dx - (radius << 1);

    while (x >= y)
    {
        color(x0 + x, y0 + y,255,255,255,img);
        color(x0 + y, y0 + x,255,255,255,img);
        color(x0 - y, y0 + x,255,255,255,img);
        color(x0 - x, y0 + y,255,255,255,img);
        color(x0 - x, y0 - y,255,255,255,img);
        color(x0 - y, y0 - x,255,255,255,img);
        color(x0 + y, y0 - x,255,255,255,img);
        color(x0 + x, y0 - y,255,255,255,img);

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
/*
void linealAngulos(StokesImage vs, int x, int y, WindowAbs * ven, CImg<float> & anterior, CImg<float> & angles){
	vector<Peso> winWeights = ven->getWeights(x,y);
	for(int auxV=0;auxV<winWeights.size();auxV++){
		Peso p=winWeights[auxV];
		float final = angles(p.f,p.c,0,0)/2*cimg::PI/180.0f;
		for(int k=0;k<3;k++){
			if(p.p>=0.05f){
				anterior(p.f,p.c,k)=p.p*(fmin(1,fmax(0,(vs.atI(p.f,p.c,k)+vs.atQ(p.f,p.c,k)*cos(2*final)+vs.atU(p.f,p.c,k)*sin(2*final))/2))*255)+(1-p.p)*anterior(p.f,p.c,k);
			}
		}
	}
}
*/

class linealAngulosPara{
	private:
		StokesImage * vs;
		bool maximizar;
		thread * vt[MAX_THREADS];
		int radio;
		vector<Peso> venPesos;
		CImg<float> angulos;
	public:
		void operator()();
		linealAngulosPara(StokesImage * vs, vector<Peso> venPesos, CImg<float> angulos){
			this->vs = vs;
			this->venPesos = venPesos;
			this->angulos=angulos;
		}

		void filtrarLineaN(int i,int n,CImg<float> & anterior){
			//cout << "PROCESO DE " << i << " A " << n << endl;
			for(int auxV=i;auxV<n;auxV++){
				if(auxV < venPesos.size()){
					Peso p=venPesos[auxV];

					if(p.p>=0.05f){
						float final = angulos(p.f,p.c,0,0)/2*cimg::PI/180.0f;
						for(int k=0;k<3;k++){
							anterior(p.f,p.c,k)=p.p*(fmin(1,fmax(0,(vs->atI(p.f,p.c,k)+vs->atQ(p.f,p.c,k)*cos(2*final)+vs->atU(p.f,p.c,k)*sin(2*final))/2))*255)+(1-p.p)*anterior(p.f,p.c,k);
						}
					}
				}
			}
		}

		thread tothread(int i, int n, CImg<float> & anterior){
			return thread(&linealAngulosPara::filtrarLineaN,this,i,n,std::ref(anterior));
		}
};

void linealAngulos(StokesImage * vs, int x, int y, WindowAbs * ven, CImg<float> & anterior, CImg<float> & angulos){
	vector<Peso> venPesos = ven->getWeights(x, y);
	linealAngulosPara fp (vs, venPesos,angulos);
	int total = 4;
	thread vt[total];
	//cout << "THREADS: " << total << " - PIXELES: " << winWeights.size()/total << endl;
	for(int i=0;i<total;i++){
		vt[i]= fp.tothread(i*venPesos.size()/total, (i+1)*venPesos.size()/total, anterior);
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
	luminanceLocalThreads fpbmax (vs,10,true);
    fpbmax.apply();
	angBmax = fpbmax.angles;

	luminanceLocalThreads fpbmin (vs,10,false);
    fpbmin.apply();
	angBmin = fpbmin.angles;
	cout << "CALCULADO ANGULOS BRILLO" << endl;
	//ANGULOS CONTRASTE
	/*
	GaussianWindow ac(vs.Stokes.width(),vs.Stokes.height(),10,5);
	WindowAbs * ptrc = &ac;
	filtroParaleloC fpcmax (vs,10,true,ptrc);
	fpcmax.filtrado();
	angCmax = fpcmax.angles;

	filtroParaleloC fpcmin (vs,10,false,ptrc);
	fpcmin.apply();
	angCmin = fpcmin.angles;
	cout << "CALCULADO ANGULOS CONTRASTE" << endl;*/
	//ANGULOS SATURACION
	GaussianWindow a2(vs.Stokes.width(),vs.Stokes.height(),15,7.5);
	WindowAbs * ptr2 = &a2;
	filtroParaleloSatVen fpsmax (vs,10,true,nullptr,ptr2,false);
	fpsmax.filtrado();
	angSmax = fpsmax.angularMean;
	cout << "CALCULADO ANGULOS SATURACION" << endl;
	filtroParaleloSatVen fpsmin (vs,10,false,nullptr,ptr2,false);
	fpsmin.filtrado();
	angSmin = fpsmin.angularMean;


	for(int i=0;i<nueva.width();i++){
		for(int j=0;j<nueva.height();j++){
			for(int k=0;k<3;k++){
				nueva(i,j,k)=vs.atI(i,j,k)/2*255;
			}
		}
	}
	CImgDisplay main_disp(nueva,"Imagen",0);
	int tam=20;
	bool maximizar=false;
	bool pulsadaM=true, newRect=true;
	Rectangulo rect;
	int modo=0; //0 brillo, 1 contraste, 2 sat, 3 original, 4 lineal
	bool cambioCon=false;
	bool pulsada = true;
	float angulo=0;
	FPS contadorfps;
	//contadorfps.empezar();
	GaussianWindow a(nueva.width(),nueva.height(),tam,tam/3.5f);
	while (!main_disp.is_closed()){
		clock_t begin = clock();
		//contadorfps.showFPS();
		CImg<float> mostrar(nueva);
		main_disp.wait();
		clock_t end = clock();
		double msi = double(end-begin)/(CLOCKS_PER_SEC/1000);
		//cout << (main_disp.button()&1) << endl;
		begin = clock();
		if(main_disp.is_keyPADADD()){
			tam++;
			cout << "INCREMENTA VENTANA" << endl;
			cout << tam << endl;
			a = GaussianWindow(nueva.width(),nueva.height(),tam,tam/3.5f);
		}
		if(main_disp.is_keyPADSUB()){
			if(tam>2){
				tam--;
				cout << "DECREMENTA VENTANA" << endl;
				a = GaussianWindow(nueva.width(),nueva.height(),tam,tam/3.5f);
			}
		}
		//MAXIMIZAR O MINIMIZAR
		if(main_disp.is_keyM()){
			if(pulsadaM){
				cambioCon=true;;
				cout << "CAMBIA MAXIMIZAR" << endl;
				maximizar = !maximizar;
				if(maximizar){
					cout << "MAXIMIZANDO" << endl;
				}
				else{
					cout << "MINIMIZAR" << endl;
				}
				pulsadaM=false;
			}
		}
		else{
			pulsadaM=true;
		}

		if(main_disp.is_keyB()){
			if(pulsada){
				cout << "OPTIMIZANDO EL BRILLO" << endl;
				modo=0;
				pulsada = false;
			}
		}
		else if(main_disp.is_keyC()){
			if(pulsada){
				cout << "OPTIMIZANDO EL CONTRASTE" << endl;
				newRect=true;
				modo=1;
				pulsada = false;
			}
		}
		else if(main_disp.is_keyS()){
			if(pulsada){
				cout << "OPTIMIZANDO LA SATURACION" << endl;
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
		else if(main_disp.is_keyF()){
			if(pulsada){
				cout << "FILTRO LINEAL" << endl;
				angulo+=10;
				modo=4;
				pulsada = false;
			}
		}
		else if(main_disp.is_keyG()){
			if(pulsada){
				cout << "Introduzca la ruta donde guardar la imagen: " << endl;
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
		if(modo==0) {
			drawcircle(main_disp.mouse_x(),main_disp.mouse_y(),tam,mostrar);
			if(main_disp.button()&1 && main_disp.mouse_y()>=0 && main_disp.mouse_x()>=0){
				//cout << "DENTRO: " << (main_disp.button()&1) << endl;
				int x = main_disp.mouse_x(), y = main_disp.mouse_y();

				//filtrosBrilloParaleloPincel(&vs, maximize, a.radio,a.getWeights(x,y), nueva);
				if(maximizar){
					linealAngulos(&vs,x,y,&a,nueva,angBmax);
				}
				else{
					linealAngulos(&vs,x,y,&a,nueva,angBmin);
				}
			}	
		}
		else if(modo==1){
			/*
			if(main_disp.button()&1 && main_disp.mouse_y()>=0 && main_disp.mouse_x()>=0){
				int x = main_disp.mouse_x(), y = main_disp.mouse_y();
				if(newRect){
					rect.i.c=x;rect.i.f=y;
					newRect=false;
				}
				rect.p.c = x; rect.p.f = y;
				rect.actu=true;
				rect.dibujar(mostrar);
			}
			if(main_disp.button()&2 && main_disp.mouse_y()>=0 && main_disp.mouse_x()>=0){
				newRect=true;
				rect.actu=false;
			}
			if(!(main_disp.button()&1) && main_disp.mouse_y()>=0 && main_disp.mouse_x()>=0){
				if(rect.actu){
					contrastePixel(vs,maximize, rect.getPixeles(), nueva);
					rect.actu=false;
				}
				newRect=true;
			}*/
			drawcircle(main_disp.mouse_x(),main_disp.mouse_y(),tam,mostrar);
			if(main_disp.button()&1 && main_disp.mouse_y()>=0 && main_disp.mouse_x()>=0){
				//cout << "DENTRO: " << (main_disp.button()&1) << endl;
				int x = main_disp.mouse_x(), y = main_disp.mouse_y();
				///cout << "FILA: " << y << " COL: " << x << endl;
				//pixeles.push_back(Pixel(y,x));
				//brilloPixel(vs,maximize,pixeles,nueva);
				//contrastePincel(&vs,maximize,x,y,&a,nueva,cambioCon);
				//filtrosConParaleloPincel(&vs,maximize, a.radio, a.getWeights(x,y), nueva, x,y,&a);

				if(maximizar){
					linealAngulos(&vs,x,y,&a,nueva,angCmax);
				}
				else{
					linealAngulos(&vs,x,y,&a,nueva,angCmin);
				}
			}
		}
		else if(modo==2){
			drawcircle(main_disp.mouse_x(),main_disp.mouse_y(),tam,mostrar);
			if(main_disp.button()&1 && main_disp.mouse_y()>=0 && main_disp.mouse_x()>=0){
				//cout << "DENTRO: " << (main_disp.button()&1) << endl;
				int x = main_disp.mouse_x(), y = main_disp.mouse_y();

				//filtrosSatParaleloPincel(&vs, maximize, a.radio,a.getWeights(x,y), nueva);
				if(maximizar){
					linealAngulos(&vs,x,y,&a,nueva,angSmax);
				}
				else{
					linealAngulos(&vs,x,y,&a,nueva,angSmin);
				}
			}
		}
		else if(modo==3){
			drawcircle(main_disp.mouse_x(),main_disp.mouse_y(),tam,mostrar);
			if(main_disp.button()&1 && main_disp.mouse_y()>=0 && main_disp.mouse_x()>=0){
				//cout << "DENTRO: " << (main_disp.button()&1) << endl;
				int x = main_disp.mouse_x(), y = main_disp.mouse_y();
				///cout << "FILA: " << y << " COL: " << x << endl;
				//pixeles.push_back(Pixel(y,x));
				//brilloPixel(vs,maximize,pixeles,nueva);
				original(vs,maximizar,x,y,&a,nueva);
			}
		}
		else if(modo==4){
			drawcircle(main_disp.mouse_x(),main_disp.mouse_y(),tam,mostrar);
			if(main_disp.button()&1 && main_disp.mouse_y()>=0 && main_disp.mouse_x()>=0){
				//cout << "DENTRO: " << (main_disp.button()&1) << endl;
				int x = main_disp.mouse_x(), y = main_disp.mouse_y();
				///cout << "FILA: " << y << " COL: " << x << endl;
				//pixeles.push_back(Pixel(y,x));
				//brilloPixel(vs,maximize,pixeles,nueva);
                linearWindow(vs, maximizar, x, y, &a, nueva, angulo);
			}
		}
		main_disp.display(mostrar);
		//contadorfps.frame();
		end = clock();
		//double msf = double(end-begin)/(CLOCKS_PER_SEC/1000);
		//cout << msi << " - " << msb << " - " << msf << endl;		
		//cout << "-------------------------" << endl;
	}
	//contadorfps.parar();
}