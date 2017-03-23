
;; * Hacked from Pipeline/Run6/define_kidpar.pro and Labtools/NP/Run6/update_kidpar_ref.pro
;;-----------------------------------------------------------------------------------------

ps       = 1
png      = 0
charsize = 1.5

;; Run6 focal planes
;; kidpar_1mm_infile = !nika.off_proc_dir+"/kidpar_ref_434pixref_1mm_bestscans_v4.fits"
;; kidpar_2mm_infile = !nika.off_proc_dir+"/kidpar_ref_434pixref_2mm_bestscans_v4.fits"

;kidpar_1mm_infile = !nika.off_proc_dir+"/kidpar_ref_nopixref_1mm_bestscans_v4.fits"
;kidpar_2mm_infile = !nika.off_proc_dir+"/kidpar_ref_nopixref_2mm_bestscans_v4.fits"

;; Get new beams for Run5
run = 5
kidpar_1mm_infile = "run5_20121123s197_1mm.fits"
kidpar_2mm_infile = "run5_20121123s197_2mm.fits"

k1  = mrdfits( kidpar_1mm_infile, 1)
k2  = mrdfits( kidpar_2mm_infile, 1)
w11 = where( k1.type eq 1, nw11)
w12 = where( k2.type eq 1, nw12)
print, nw11, nw12

fwhm_1mm = sqrt( k1[w11].fwhm_x*k1[w11].fwhm_y)
fwhm_2mm = sqrt( k2[w12].fwhm_x*k2[w12].fwhm_y)
str = {a:fwhm_1mm, b:fwhm_2mm}

binsize = 0.2
blend = 1

wind, 1, 1, /free, /large
np_histo, str, xhist, yhist, gpar, blend=blend, fcol=[70,250], /nolegend, $
          xtitle='Beam FWHM [Arcsec]', /fit, xrange=[10, 22], binsize=binsize, yrange=[0,30], $
          charsize=charsize
if (png eq 1) or (ps eq 1) then outplot, file='run5_beam_stats', /png, charsize=1.5, thick=2, charthick=2
np_histo, str, xhist, yhist, blend=blend, fcol=[70,250], /nolegend, $
          xtitle='Beam FWHM [Arcsec]', xrange=[10, 22], binsize=binsize, yrange=[0,30], $
          charsize=charsize
legendastro, ['Run 5', '', $
              'Average FWHM = '+num2string(gpar[0,1])+" +- "+num2string(gpar[0,2]), '', $
              'Average FWHM = '+num2string(gpar[1,1])+" +- "+num2string(gpar[1,2])], $
             box=0, chars=charsize, textcol=[0,0,70,0,250]
;legendastro, ['', '', ], box=0, chars=charsize, textcol=250, /right
outplot, /close

if (ps eq 1) then spawn, 'convert run5_beam_stats.png run5_beam_stats.eps'
;;image = TVREAD(/true)
;;image = transpose(image,[1,2,0])
;;red = image[*,*, 0]
;;grn = image[*,*, 1]
;;blu = image[*,*, 2]
;;whiteIndices = Where((red eq 255) and (grn eq 255) and (blu eq 255), count)
;;s = Size(image, /DIMENSIONS)
;;xsize = s[0] & ysize = s[1]
;;alpha = BytArr(xsize, ysize) + 255B
;;IF count GT 0 THEN alpha[whiteIndices] = 0
;;transparentImage1 = [ [[red]], [[grn]], [[blu]], [[alpha]] ]
;;alpha1 = alpha
;;red1 = red
;;blu1 = blu
;;grn1 = grn

;; Get new beams for Run6
kidpar_1mm_infile = "my_new_ref_run6_1mm.fits"
kidpar_2mm_infile = "my_new_ref_run6_2mm.fits"

k1  = mrdfits( kidpar_1mm_infile, 1)
k2  = mrdfits( kidpar_2mm_infile, 1)
w11 = where( k1.type eq 1, nw11)
w12 = where( k2.type eq 1, nw12)
print, nw11, nw12

fwhm_1mm = sqrt( k1[w11].fwhm_x*k1[w11].fwhm_y)
fwhm_2mm = sqrt( k2[w12].fwhm_x*k2[w12].fwhm_y)
str = {a:fwhm_1mm, b:fwhm_2mm}
wind, 1, 1, /free, /large
np_histo, str, xhist, yhist, gpar, blend=blend, fcol=[70,250], /nolegend, $
          xtitle='Beam FWHM [Arcsec]', /fit, xrange=[10, 22], binsize=binsize, yrange=[0,30], $
          charsize=charsize
if (ps eq 1) or (png eq 1) then outplot, file='run6_beam_stats', /png, charsize=1.5, thick=2, charthick=2
np_histo, str, xhist, yhist, blend=blend, fcol=[70,250], /nolegend, $
          xtitle='Beam FWHM [Arcsec]', xrange=[10, 22], binsize=binsize, yrange=[0,30], $
          charsize=charsize
legendastro, ['Run 6', '', $
              'Average FWHM = '+num2string(gpar[0,1])+" +- "+num2string(gpar[0,2]), '',$
              'Average FWHM = '+num2string(gpar[1,1])+" +- "+num2string(gpar[1,2])], $
             box=0, chars=charsize, textcol=[0,0,70,0,250]
outplot, /close
if (ps eq 1) then spawn, 'convert run6_beam_stats.png run6_beam_stats.eps'

;if ps eq 1 or png eq 1 then spawn, 'convert run'+strtrim(run,2)+'_beam_stats.png run'+strtrim(run,2)+'_beam_stats.eps'
;if ps eq 1 or png eq 1 then spawn, 'open run'+strtrim(run,2)+'_beam_stats.eps &'
;; image = TVREAD(/true)
;; image = transpose(image,[1,2,0])
;; red = image[*,*, 0]
;; grn = image[*,*, 1]
;; blu = image[*,*, 2]
;; whiteIndices = Where((red eq 255) and (grn eq 255) and (blu eq 255), count)
;; s = Size(image, /DIMENSIONS)
;; xsize = s[0] & ysize = s[1]
;; alpha = BytArr(xsize, ysize) + 255B
;; IF count GT 0 THEN alpha[whiteIndices] = 0
;; transparentImage2 = [ [[red]], [[grn]], [[blu]], [[alpha]] ]
;; alpha2 = alpha
;; red2 = red
;; blu2 = blu
;; grn2 = grn
;; 
;; alpha1 = double(alpha1)
;; alpha2 = double(alpha2)
;; norm   = double(alpha1) + double(alpha2)
;; w = where( norm ne 0)
;; alpha1[w] = alpha1[w]/norm[w]
;; alpha2[w] = alpha2[w]/norm[w]
;; 
;; red_merge = alpha1*red1 + alpha2*red2
;; blu_merge = alpha1*blu1 + alpha2*blu2
;; grn_merge = alpha1*grn1 + alpha2*grn2
;; 
;; red_merge = byte(red_merge/max(red_merge) * 255)
;; blu_merge = byte(blu_merge/max(blu_merge) * 255)
;; grn_merge = byte(grn_merge/max(grn_merge) * 255)
;; alpha_merge = byte( (norm gt 0)*255)
;; image_merge = [ [[red_merge]], [[grn_merge]], [[blu_merge]], [[alpha_merge]] ]
;; 
;; transparentImage = Transpose(image_merge, [2,0,1])
;; wind, 3, 3, /free
;; cgimage, transparentimage

end
