#include <iostream>
#include <chrono>
#include "Stokes/Stokes.h"
#include "Stokes/LuminanceFilters.h"
#include "Stokes/ContrastFilters.h"
#include "Stokes/SaturationFilters.h"
#include "Stokes/Window.h"

#include <ctime>
#include <string>

using namespace cimg_library;
using namespace std;



int main(int argc, char* argv[]) {
	string dir="";
	string out="";
	int grados=0,delta=0;
	int maximizar=-5;
	bool global=false;
	int tipoFiltro = -1;
	bool ventana=false;
	int radio=0;
	float desviacion=0;
	bool saveVS=false;
	string fichVS="";
	bool rcp = false, lcp=false, param = false;
	float c1, c2, c3;

	bool blur = false;
	float desvB = 0;
	int radioB=0;

	bool circular=false;

	StokesImage * vs;
	if(argc<3){
		std::cerr << "Usage: " << argv[0] << " -i stokesDirectory -o imgResul " << endl;

			std::cerr << "[-linear degrees] Linear polarizer" << endl;
			std::cerr << "[-rcp] Right circular polarizer" << endl;
			std::cerr << "[-lcp] Left circular polarizer" << endl;
			std::cerr << "[-circular alpha delta] " << endl;
			std::cerr << "[-param c1 c2 c3] Custom polarizer" << endl;
			std::cerr << "[-l (minimize|maximize) 0|1] Local Luminance " << endl;
			std::cerr << "[-lg (minimize|maximize) 0|1] Global Luminance " << endl;
			std::cerr << "[-cv (minimize|maximize) 0|1 radius standardDeviation] Local Contrast " << endl;
			std::cerr << "[-cg (minimize|maximize) 0|1] Global Contrast " << endl;
			std::cerr << "[-s (minimize|maximize) 0|1] Local Saturation" << endl;
			std::cerr << "[-sv (minimize|maximize) 0|1 radius standardDeviation] Local Saturation with window" << endl;
			std::cerr << "with -s and -sv, [-blur radius standardDeviation]to blur the angle" << endl;
			std::cerr << "[-sg (minimize|maximize) 0|1] Global Saturation " << endl;
			std::cerr << "for more info read the readme.txt" << endl;
		return -1;
	}
	for(int i=0;i<argc;i++){
		string arg = argv[i];
		//DIRECTORIO DE ENTRADA
		if(arg=="-i"){
			i++;
			dir=argv[i];
			clock_t begin = clock();
			vs = new StokesImage(dir);
			clock_t end = clock();
		 	double msec = double(end - begin) / (CLOCKS_PER_SEC/1000);
		 	cout << "StokesImage created in: " << msec << endl;
		}
		//IMAGEN DE STOKES DE ENTRADA
		else if(arg=="-vsi"){
			i++;
			string auxVS=argv[i];
			vs = new StokesImage();
			vs->load(auxVS);
		}

		//GUARDAR LOS VECTORES EN UNA IMAGEN .cimg
		else if(arg=="-vso"){
			saveVS = true;
			i++;
			fichVS = argv[i];
		}
		//IMAGEN DE SALIDA
		else if(arg=="-o"){
			i++;
			out=argv[i];
		}
		//GRADOS
		else if(arg=="-lineal"){
			i++;
			maximizar=0;
			string g=argv[i];
			grados=stoi(g);
		}
		//RCP
		else if(arg=="-rcp"){
			maximizar=0;
			rcp = true;
		}
		//LCP
		else if(arg=="-lcp"){
			maximizar=0;
			lcp = true;
		}
		else if(arg=="-param"){
			
			maximizar=0;
			param=true;
			i++;
			string g=argv[i];
			c1=stoi(g);
			i++;
			g=argv[i];
			c2=stoi(g);
			i++;
			g=argv[i];
			c3=stoi(g);
		}
		//LCP
		else if(arg=="-circular"){
			maximizar=0;
			circular=true;
			i++;
			string g=argv[i];
			grados=stoi(g);
			i++;
			g=argv[i];
			delta=stoi(g);
		}
		//BRILLO PIXEL
		else if(arg=="-l"){
			tipoFiltro=0;
			i++;
			string maxi=argv[i];
			if(maxi=="0"){
				maximizar=-1;
			}
			else{
				maximizar=1;
			}
		}
		//BRILLO GLOBAL
		else if(arg=="-lg"){
			tipoFiltro=0;
			global=true;
			i++;
			string maxi=argv[i];
			if(maxi=="0"){
				maximizar=-1;
			}
			else{
				maximizar=1;
			}
		}
		//CONTRASTE PIXEL
		else if(arg=="-c"){
			tipoFiltro=1;
			i++;
			string maxi=argv[i];
			if(maxi=="0"){
				maximizar=-1;
			}
			else{
				maximizar=1;
			}
		}
		//CONTRASTE GLOBAL
		else if(arg=="-cg"){
			tipoFiltro=1;
			global=true;
			i++;
			string maxi=argv[i];
			if(maxi=="0"){
				maximizar=-1;
			}
			else{
				maximizar=1;
			}
		}
		//CONTRASTE VENTANA
		else if(arg=="-cv"){
			tipoFiltro=1;
			global=false;
			ventana=true;
			i++;
			string maxi=argv[i];
			if(maxi=="0"){
				maximizar=-1;
			}
			else{
				maximizar=1;
			}
			i++;
			string radioS=argv[i];
			radio=stoi(radioS);
			i++;
			string desvS=argv[i];
			desviacion=stof(desvS);
		}
		//SATURACION PIXEL
		else if(arg=="-s"){
			tipoFiltro=2;
			i++;
			string maxi=argv[i];
			if(maxi=="0"){
				maximizar=-1;
			}
			else{
				maximizar=1;
			}
		}
		//SATURACION GLOBAL
		else if(arg=="-sg"){
			tipoFiltro=2;
			global=true;
			i++;
			string maxi=argv[i];
			if(maxi=="0"){
				maximizar=-1;
			}
			else{
				maximizar=1;
			}
		}
		else if(arg=="-sv"){
			tipoFiltro=2;
			global=false;
			ventana=true;
			i++;
			string maxi=argv[i];
			if(maxi=="0"){
				maximizar=-1;
			}
			else{
				maximizar=1;
			}
			i++;
			string radioS=argv[i];
			radio=stoi(radioS);
			i++;
			string desvS=argv[i];
			desviacion=stof(desvS);
		}
		else if(arg=="-blur"){
			blur=true;
			i++;
			string radioS=argv[i];
			radioB=stoi(radioS);
			i++;
			string desvS=argv[i];
			desvB=stof(desvS);
		}
	}
	if(saveVS){
		vs->save(fichVS);
	}
	clock_t begin,end;
	double msec;

	//StokesImage vsP(CImg<float>("prueba8/salida/I.jpg"),CImg<float>("prueba8/salida/Q.jpg"),CImg<float>("prueba8/salida/U.jpg"),CImg<float>("prueba8/salida/V.jpg"));
	//vs = new StokesImage("escenas/I.hdr","escenas/Q.hdr","escenas/U.hdr","escenas/V.hdr");
	/*CImgDisplay S0(vs.I, "I");
	S0.resize(720,540);
	CImgDisplay S1(vs.Q, "Q");
	S1.resize(720,540);
	CImgDisplay S2(vs.U, "U");
	S2.resize(720,540);
	CImgDisplay S3(vs.V, "V");
	S3.resize(720,540);*/
	//CImg<float> filtrada=filtroLinealPorPixel(vs,false);
	/* nueva(filtrada, "xPIXEL");
	nueva.resize(720,540);*/
	vs->resize(vs->Stokes.width()*1.0f/vs->Stokes.height()*720,720);
	if(maximizar==0){
		begin = clock();
		CImg<float> filtrolineal;
		if(rcp){
			filtrolineal= RCP(*vs)*255;
		}
		else if(lcp){
			filtrolineal= LCP(*vs)*255;
		}
		else if(param){
			filtrolineal= customFilter(*vs, c1, c2, c3)*255;
		}
		else if(circular){
			filtrolineal= circularFilter(*vs, deg2rad(grados), deg2rad(delta))*255;
		}
		else{
			filtrolineal= linearFilter(*vs, grados)*255;
		}
		end = clock();
		msec = double(end - begin) / (CLOCKS_PER_SEC/1000);
		cout << "Execution time: " << msec << endl;
		filtrolineal.save_jpeg(out.c_str());
	}
	else if(tipoFiltro == 0){
		bool max= maximizar==1;
		begin = clock();
		CImg<float> filtrolineal;
		if(!global){
			//std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

			filtrolineal= localLuminance(*vs, max)*255;

			//std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();

			//std::cout << "Time difference PARALELO = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms" <<std::endl;

		}
		else{
			filtrolineal= globalLuminance(*vs, max)*255;
		}
		end = clock();
		msec = double(end - begin) / (CLOCKS_PER_SEC/1000);
		cout << "ORIGINAL LUMINANCE: " << luminanceIndex(vs->calculateI() / 2 * 255) << endl;
		cout << "FILTERED LUMINANCE: " << luminanceIndex(filtrolineal) << endl;

		cout << "Execution time: " << msec << "ms" << endl;
		filtrolineal.save_jpeg(out.c_str());
	}
	else if(tipoFiltro == 1){
		bool max= maximizar==1;
		begin = clock();
		CImg<float> filtrolineal;
		if(!ventana){
			if(!global){
				BasicWindow a(vs->Stokes.width(),vs->Stokes.height());
				filtrolineal= localContrast(*vs, max, &a)*255;
			}
			else{
				filtrolineal= globalContrast(*vs,max)*255;
			}
		}
		else{
			GaussianWindow a(vs->Stokes.width(),vs->Stokes.height(),radio,desviacion);
			cout << "WINDOW RADIUS = " << radio << " STANDARD DEVIATION = " << desviacion << endl; 
			WindowAbs * ptr = &a;
			filtrolineal=localContrast(*vs,max,ptr)*255;
		}

		cout << "ORIGINAL CONTRAST: " << contrastIndex(vs->calculateI()/2) << endl;
		cout << "FILTERED CONTRAST: " << contrastIndex(filtrolineal*1.0f/255) << endl;
		end = clock();
		msec = double(end - begin) / (CLOCKS_PER_SEC/1000);
		cout << "Execution time: " << msec << "ms" << endl;
		filtrolineal.save_jpeg(out.c_str());
	}
	else if(tipoFiltro == 2){
		bool max= maximizar==1;
		begin = clock();
		CImg<float> filtrolineal;
		if(!ventana){
			if(!global){
				WindowAbs * ptr2 = nullptr;
				GaussianWindow b(vs->Stokes.width(),vs->Stokes.height(),radioB,desvB);
				if(blur){
					ptr2=&b;
				}
				filtrolineal=localSaturation(*vs,max,nullptr,ptr2)*255;
			}
			else{
				filtrolineal=globalSaturation(*vs,max)*255;
			}
		}
		else{
			GaussianWindow a(vs->Stokes.width(),vs->Stokes.height(),radio,desviacion);
			cout << "WINDOW RADIUS = " << radio << " STANDARD DEVIATION = " << desviacion << endl; 
			WindowAbs * ptr = &a;
			//filtrolineal=saturacionVentana(*vs,max,ptr)*255;
			WindowAbs * ptr2 = nullptr;
			GaussianWindow b(vs->Stokes.width(),vs->Stokes.height(),radioB,desvB);
			if(blur){
				ptr2=&b;
			}
			filtrolineal=localSaturation(*vs,max,ptr,ptr2)*255;
		}
		cout << "ORIGINAL SATURATION: " << saturationIndex(vs->calculateI()/2*255) << endl;
		cout << "FILTERED SATURATION: " << saturationIndex(filtrolineal) << endl;
		end = clock();
		msec = double(end - begin) / (CLOCKS_PER_SEC/1000);
		cout << "Execution time: " << msec << "ms" << endl;
		filtrolineal.save_jpeg(out.c_str());
	}

	/*begin = clock();
	CImg<float> filtrolinealP=filtroLinealC(vsP,135);
	end = clock();
	msec = double(end - begin) / (CLOCKS_PER_SEC/1000);
	cout << "FiltroLinealMALO en: " << msec << endl;*/
	/*CImgDisplay lineal(filtrolineal, "LINEAL");
	lineal.resize(720,540);*/
	/*
	CImg<float> asaco=vs.filtroASaco(false);
	CImgDisplay asacoA(asaco, "FILTRO A SACO");
	asacoA.resize(720,540);
	vs.I.save("I.jpg");
	vs.Q.save("Q.jpg");
	vs.U.save("U.jpg");
	vs.V.save("V.jpg");*/
	//filtrada.normalize(0,255);
	
	//filtrolinealP.normalize(0,255);
	//filtrada.save("xPIXELnoReflejoNUEVO.jpg");
	//filtrolinealP.save("linealCHUNGO.jpg");
	/*asaco.save("parametrosnoReflejo.jpg");*/
	return 0;
}