; Put into the NikaDataProductsVB.0.tex

The IRAM NIKA OpenPool2 products v0
------------------------------------

This version V0 of the NIKA second open pool products has been made by the NIKA
team by the 5 December 2014. 
It is distributed by IRAM. It follows the
observing run by 2 weeks which were needed for the data processing.

The products are described in this document:

http://www.iram.es/IRAMES/mainWiki/Continuum/NIKA/DataReduction?action=AttachFile&do=view&target=NikaDataProductsVB.0.pdf

For each project, there are fits files containing the maps of each scan and
combined scans of each object. The data processing has been oriented towards
the diffuse emission. 

There are also illustrative pdf files. You will find that for each object,
there is a directory that contains:

- The fits files (RaDec2000 projection) for target sources

- Figures: many pdf files were generated out of the fits files to get a quick
look at them. Units in Jy/beam unless otherwise stated. A smoothing has been
done with a Gaussian FWHM of 10 arcseconds (only for the figures). 

JK pdf are jackknife maps representing the half difference of a random split
of scans in two halves.

SNR pdf is the signal-to-noise maps assuming Gaussian white noise.

flux pdf is the brightness map display

stddev pdf is the standard deviation map display

time pdf is the total integration time per 2 arcsecond square pixel.

scan pdf is a brightness display per scan


The main beam calibration was done assuming a Gaussian main beam of 12.5 and
18.5 arcseconds (FWHM). The primary calibrator is Uranus with fluxes of 43.0
and 18.0 Jy.

The main beam to full beam correction is by
1.56+-0.10 at 1mm and 1.35+-0.10 at 2mm. It has not been applied to the maps.

At this stage, the offline products cannot be used for scientific analysis as
some photometric uncertainties are still there (mostly opacity effects). The
point-source photometry may probably correct at the 20% level at 1mm and 15%
level at 2mm. These products must just be used to evaluate the potential
return of each observed source.

We hope to deliver a final V1 version of the second NIKA open pool data
products in several months including our best strategy of systematic removals.
 
Contact your NIKA friend of project to give us your feedback. 

The NIKA team 
