NIKA products v1
----------------

This version V1 of the NIKA first open pool products has been made by the NIKA
team by the end of May 2014. It is distributed by IRAM. It follows the quick
release done at the end of the pool (V0).

The products are described in this document:

http://www.iram.es/IRAMES/mainWiki/Continuum/NIKA/DataReduction?action=AttachFile&do=view&target=NikaDataProductsV1.pdf

For each project, there are fits files containing the maps of each scan and
combined scans of each object. The data processing has been oriented towards
point-source extraction or diffuse emission depending on the aim of your
project.  The main processing has improved by including more scans and dealing
with internal problems more effectively.

There are also logbook files and illustrative pdf files. On the main directory
you will find:

- logfile_all_scans.csv : is a comma separated file ready for your favorite
spreadsheet handler (libreoffice or excel) It gives the summary information of
all the scans that could be used including opacity measurements, sky noise in
Jy (per band of frequency from 1mHz to 10 Hz), camera noise NEFD (in
mJy.s^1/2) and finally a flux in mJy of a point-source at the center of the
map (fixed central position, fixed FWHM, Gaussian fit) with the error at both
wavelengths. For pointing sources the position is allowed to vary.

Then for each object, there is a directory that contains:

- The fits files (RaDec2000 projection for target sources, AzEl projection for
  pointing sources and Planets)

- flux_source.csv : a comma separated file giving the flux as measured in the
scan maps and the combined scans maps.

- Figures: many pdf files were generated out of the fits files to get a quick
look at them. Units in Jy/beam unless otherwise stated.

The calibration was done assuming a Gaussian main beam of 12.5 and 18.5
arcseconds (FWHM). The primary calibrator is Uranus with fluxes of 36.4 and
15.3 Jy.

For Lissajous scans (V1Liss directories), the pipeline introduces a filtering
effect for point sources of about 0.70 and 0.56 (at 1 and 2mm). This effect
has been corrected in the .csv files but not in the maps. Thus fluxes in Jy in
the .csv files should be correct.  

For you to assess the filtering effect, we have rerun the pipeline with a much
less agressive decorrelation that should not affect the fluxes. For that just
look at V1cLiss directories made for calibration sources only.

For OTF maps (directories V1otf): the main beam to full beam correction is by
1.56+-0.10 at 1mm and 1.35+-0.10 at 2mm. It has not been applied to the maps.

At this stage, the offline products can be used for a first scientific
analysis as most of the scans have been included. The point-source photometry
is probably correct at the 15% level at 1mm and 10% level at 2mm.

We hope to deliver a final V2 version of the first NIKA open pool data
products in the Summer including our best strategy of systematic removals.
 
Contact
your NIKA friend of project to give us your feedback. 

The NIKA team 
