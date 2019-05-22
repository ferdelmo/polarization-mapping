#ifndef STOKES_H_
#define STOKES_H_

#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include "../CImg/CImg.h"

using namespace cimg_library;
using namespace std;


/*
    Check if a files exists
 */
inline bool exists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

/*
    Class to store a Stokes image
 */
class StokesImage{
	public:
		CImg<float> Stokes;
		StokesImage(){}
		/*
			Genearte the stokes image. Recives as input a path to a directory with the 5 necesary images:
		    im0.jpg, im45.jpg, im90.jpg, im135.jpg and imcd.jpg (images with a linear polarizer at 0,45,
		    90, 135 degrees and the image with the circular polarizer
		*/
		StokesImage(string dir) {
            CImg<float> im0;
            im0.load_jpeg((dir + "/im0.jpg").c_str());
            CImg<float> im90;
            im90.load_jpeg((dir + "/im90.jpg").c_str());
            CImg<float> im45;
            im45.load_jpeg((dir + "/im45.jpg").c_str());
            CImg<float> im135;
            im135.load_jpeg((dir + "/im135.jpg").c_str());
            CImg<float> imd(im0.width(), im0.height(), 1, 3, 0);
            if (exists((dir + "/imcd.jpg").c_str())) {
                imd.load_jpeg((dir + "/imcd.jpg").c_str());
            }
            Stokes = CImg<float>(im0.width(), im0.height(), 4, 3, 0);
            //cout << im90.width() << " " << im90.height() << " " << im90.depth() << " " << im90.spectrum() << endl;
            for (int i = 0; i < im0.width(); i++) {
                for (int j = 0; j < im0.height(); j++) {
                    for (int k = 0; k < 3; k++) {
                        Stokes(i,j,0,k) = (im0(i,j,k) / 255.0f + im90(i,j,k) / 255.0f + im45(i,j,k) / 255.0f + im135(i,j,k) / 255.0f)/2;
                        Stokes(i,j,1,k) = im0(i,j,k) / 255.0f - im90(i,j,k) / 255.0f;
                        Stokes(i,j,2,k) = im45(i,j,k) / 255.0f - im135(i,j,k) / 255.0f;
                        Stokes(i,j,3,k) = 2*imd(i,j,k) / 255.0f - Stokes(i,j,0,k);

                    }
                }
            }
        }

		/*
			Aplica un linearFilter de grados g
		*/

		float atI(int x,int y, int k){
		    return Stokes(x,y,0,k);
		}
		float atQ(int x,int y, int k){
		    return Stokes(x,y,1,k);
		}
		float atU(int x,int y, int k){
		    return Stokes(x,y,2,k);
		}
		float atV(int x,int y, int k){
		    return Stokes(x,y,3,k);
		}
		/*
		    Return the parameter I as an image
		 */
		CImg<float> calculateI(){
			CImg<float> resul(Stokes.width(),Stokes.height(),1,3,0);
			for(int i=0;i<resul.width();i++){
				for(int j=0;j<resul.height();j++){
					for(int k=0;k<3;k++){
						resul(i,j,k)=atI(i,j,k);
						if(atI(i,j,k)>1){
							//cout << "Me suisido" << " " << atI(i,j,k) << endl;
						}
					}
				}
			}
			return resul;
		}
		/*
		    Return the parameter Q as an image
		 */
		CImg<float> calculateQ(){
			CImg<float> resul(Stokes.width(),Stokes.height(),1,3,0);
			for(int i=0;i<resul.width();i++){
				for(int j=0;j<resul.height();j++){
					for(int k=0;k<3;k++){
						resul(i,j,k)=(atQ(i,j,k)+1)/2;
						//resul(i,j,k)=fmax(0,(atQ(i,j,k)/atI(i,j,k))+1);
						if(resul(i,j,k)>1){
							//cout << "Me suisido" << " " << resul(i,j,k) << endl;
						}
					}
				}
			}
			return resul;
		}
		/*
		    Return the parameter U as an image
		 */
		CImg<float> calculateU(){
			CImg<float> resul(Stokes.width(),Stokes.height(),1,3,0);
			for(int i=0;i<resul.width();i++){
				for(int j=0;j<resul.height();j++){
					for(int k=0;k<3;k++){
						resul(i,j,k)=(atU(i,j,k)+1)/2;
						//resul(i,j,k)=fmax(0,(atU(i,j,k)/atI(i,j,k))+1);
						
					}
				}
			}
			return resul;
		}
		/*
		    Return the parameter V as an image
		 */
		CImg<float> calculateV(){
			CImg<float> resul(Stokes.width(),Stokes.height(),1,3,0);
			for(int i=0;i<resul.width();i++){
				for(int j=0;j<resul.height();j++){
					for(int k=0;k<3;k++){
						resul(i,j,k)=(atV(i,j,k)+1)/2;
						//resul(i,j,k)=fmax(0,(atV(i,j,k)/atI(i,j,k)+1));
						if(resul(i,j,k)>1 || resul(i,j,k)<0){
							cout << resul(i,j,k) << endl;
							resul(i,j,k)=1;
						}
					}
				}
			}
			return resul;
		}

		void resize(int width, int height){
			Stokes.resize(width,height);
		}

		void save(string f){
			Stokes.save_cimg(f.c_str());
		}
		void load(string f){
			Stokes = Stokes.load_cimg(f.c_str());
		}
};

#endif