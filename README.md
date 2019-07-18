# Polarization Mapping
A set of tools for editing polarized images and applying polarizing filters by software, which range from any commercial linear or circular filter to different optimization filters, which maximize or minimize different values of the image (luminance, contrast or saturation). These optimization filters are subdivided in two kinds:
- Global: select the best linear polarizer to maximize or minimize the parameter in the whole image.
- Local: select a different linear polarizer per pixel to maximize (or minimize) the parameter in that pixel.

Two different tools are provided:
- A command line tool, that given a Stokes image, apply the selected filter to a result image.
- An interactive tool, that allow to edit every pixel of the image separately, being able to apply both commercial and local optimization filters in a specific area using a brush.

The tool recive as input a directoy with 5 images with different polarizations filters. 4 with a linear polarizer at 0º, 45º, 90º and 135º, and another with a RCP(Right Circular Polarizer). For the optimization filters or the linear filters, the circular polarizer is not needed and the tool can work with only the 4 images with linear polarization.

We provide a [dataset](http://dx.doi.org/10.17632/s4sc8zt4sx.1) with the 5 images needed for different scenes (some scenes do not have the circular image): 

## Compile
- git: for downloading the source code
- A C++11 capable compiler.
- [CImg](https://github.com/dtschump/CImg): A open-source C++ toolkit for image processing
- [Quartic](https://github.com/sasamil/Quartic): A library to resolve equation of 4th order

```
git clone https://github.com/ferdelmo/polarization-mapping.git
cd polarization-mapping
git clone https://github.com/dtschump/CImg
git clone https://github.com/sasamil/Quartic
```
To compile in linux:
```
g++ -std=c++11  -o pm main.cc Quartic/quartic.cpp -O2 -L/usr/X11R6/lib -lm -lpthread -lX11
g++ -std=c++11  -o BrushTool BrushTool.cc Quartic/quartic.cpp -O2 -L/usr/X11R6/lib -lm -lpthread -lX11
```
Windows(MingW windows version):
```
g++ -std=c++11 -o pm.exe main.cc Quartic/quartic.cpp -O2 -lgdi32
g++ -std=c++11 -o BrushTool.exe BrushTool.cc Quartic/quartic.cpp -O2 -lgdi32
```
Mac OS X:
```
g++ -o pm.exe main.cc Quartic/quartic.cpp -O2 -lm -lpthread -I/usr/X11R6/include -L/usr/X11R6/lib -lm -lpthread -lX11
g++ -o BrushTool.exe BrushTool.cc Quartic/quartic.cpp -O2 -lm -lpthread -I/usr/X11R6/include -L/usr/X11R6/lib -lm -lpthread -lX11
```

Execution instructions can be found in the readme.txt file.