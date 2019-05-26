# Polarization Mapping
A set of tools for editing polarized images and applying polarizing filters by software, which range from any commercial linear or circular filter to different optimization filters, which maximize or minimize different values of the image (luminance, contrast or saturation). These optimization filters are subdivided in two kinds:
- Global: select the best linear polarizer to maximize or minimize the parameter in the whole image.
- Local: select a different linear polarizer per pixel to maximize (or minimize) the parameter in that pixel.

Two different tools are provided:
- A command line tool, that given a Stokes image, apply the selected filter to a result image.
- An interactive tool, that allow to edit every pixel of the image separately, being able to apply both commercial and local optimization filters in a specific area using a brush.
