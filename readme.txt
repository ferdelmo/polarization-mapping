POLARIZATION MAPPING

Compile:

g++ -std=c++11  -o pm main.cc Quartic/quartic.cpp-O2 -L/usr/X11R6/lib -lm -lpthread -lX11
g++ -std=c++11  -o BrushTool BrushTool.cc Quartic/quartic.cpp-O2 -L/usr/X11R6/lib -lm -lpthread 

Use:

We are going to use the following abbreviations:
pt -> path to the directory with the images to generate the Stokes image. This images have to be
named im0.jpg, im45.jpg, im90.jpg, im135.jpg for the images with linear polarization at 0º, 45º,
90º and 135º respectively. The circular image, if exists, has to be named imcd.jpg.

resul.jpg -> Output image

r -> radius in pixel of a window

d -> standard deviation to weight each pixel as a function of its distance to the center and a
 gaussian distribution with that standard deviation

Basic filters:

To apply a linear polarizer at x degrees:
./pm -i pt -o resul.jpg -linear x

To apply a RCP and a LCP (Right/Left circular polarizer):
./pm -i pt -o resul.jpg -rcp
./pm -i pt -o resul.jpg -lcp

To apply a elipitical polarizer composed by a linear polarizer at x degrees and a retarder of y
 degrees:
./pm -i pt -o resul.jpg -circular x y

To apply a custom parameter filter, with coefficients c1 for Q, c2 for U and c3 for V (coefficient must be between -1 and 1, and the resulting image will be I=I+c1*Q+c2*U+c3*V):
./pm -i pt -o resul.jpg -param c1 c2 c3

Global optimization filters:

To maximize the luminance:
./pm -i pt -o resul.jpg -lg 1
To minimize the luminance:
./pm -i pt -o resul.jpg -lg 0

To maximize the contrast:
./pm -i pt -o resul.jpg -cg 1
To minimize the contrast:
./pm -i pt -o resul.jpg -cg 0

To maximize the saturation:
./pm -i pt -o resul.jpg  -sg 1
To minimize the saturation:
./pm -i pt -o resul.jpg  -sg 0

Local optimization filters:
All this filters will generate the image "angles.jpg", where the color of each pixel correspond
to the angle of the linear polarizer used in that pixel, using the HSV color space.

To maximize the luminance:
./pm -i pt -o resul.jpg -l 1
To minimize the luminance:
./pm -i pt -o resul.jpg -l 0

To maximize the contrast:
./pm -i pt -o resul.jpg -c 1 r d
To minimize the contrast: 
./pm -i pt -o resul.jpg -c 0 r d
r and d define the gaussian window used to calculate the contrast in each pixel.
The radius of the window is a critical parameter, with sizes higher than ~15 the excution can be pretty slow.

To maximize the saturation:
./pm -i pt -o resul.jpg  -s 1
To minimize the saturation:
./pm -i pt -o resul.jpg  -s 0
This filter apply a different filter to each pixel, but the result can be noisy due to the variance in the angle selection of nearby pixels.
A blur in the angle selected using a gaussian window can be applied as:
./pm -i pt -o resul.jpg  -s 0 -blur r d
Again, the radius is a critical parameter, with sizes higher than ~15-20 the exuction can be slow.

Brush tool:
./BrushTool -i pt 
This command will open a window that will show the image without any polarizer and a brush.
This tool is still a research prototype and has no graphical interface. 
When you launch it, you will have to wait a bit, as the angles of the local contrast and saturation optimizacion has to be precalculated (it cannot be computed in real time).
To change the brush, enlarge or shrink it, etc. you have to press different keys:
- Press '+' and '-' to adjust the brush size
- 'l' to use the local luminance brush
- 'c' to use the local contrast brush
- 's' to use the local saturation brush
- 'm' to toggle between maximization and minimization, in maximization mode the brush will be red,
	while in minimization mode the brush will be blue.
- 'r' to use a RCP filter
- 'l' to use a LCP filter
- 'g' to save the actual image: you will have to use the console where you launched the tool to enter the path to save the image: ..../resul.jpg and press enter
- 'p' to use a linear polarizer brush: again, will have to use the console where you launched the tool to enter the degrees of the linear polarizer: 15 and press enter
- 'o' work as a polarization eraser, this brush restores the original image (the image without
 any polarizer)
