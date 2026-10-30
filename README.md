# plot_fits
Tool to obtain and plot a position-velocity diagram from a cube under a two-dimensional CASA region.

Packages needed: astropy, numpy, scipy, matplotlib, Regions (https://pypi.org/project/regions/).

This code was employed to create the plots in Martínez-Henares et al. (2025) https://ui.adsabs.harvard.edu/abs/2025A%26A...699A.382M/abstract . If you use this code for your research, a mention in your publication is much appreciated.

An example of the type of images that can be created is shown below.

![]()

## How to use the code

To use this code, first the options have to be edited in the code, and then run. Typically you should only be editing lines 267 and below.

The `area_pv.py` file is the code itself, and it includes a working example of how it should be used. The best way to make sure that you can use the code is to reuse the example I provide with your own files and custom options. 

The workflow to use this code is:

1. Introduce the path to the FITS and region files.
2. Define the center of your PV diagram.
3. Define the systemic velocity of your source.
4. Optional: decide the range of channels in the cube that will be used (i.e. upon previous check with CARTA).
5. Optional: edit the colormap. Default is a custom version of the 'turbo' matplotlib map.

When all options are set, you can run the code by

```
python area_pv.py
```

## Potential issues

Very large areas will not produce meaningful PV diagrams. 
This code was thought for outflows. Therefore, it has been tested mostly with elongated-shaped areas.
This code was developed for a very specific case. To solve potential issues, some level of user editing would be necessary. 

## Bugs, questions

Feel free to contact me at antonio.martinezhenares@riken.jp . 

If you use this code for your research and you produce figures in your paper, I would greatly appreciate a citation on Martínez-Henares et al. (2025) https://ui.adsabs.harvard.edu/abs/2025A%26A...699A.382M/abstract .
