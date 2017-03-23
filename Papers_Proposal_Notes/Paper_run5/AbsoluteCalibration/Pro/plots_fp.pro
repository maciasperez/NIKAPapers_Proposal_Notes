
;; * Hacked from Pipeline/Run6/define_kidpar.pro and Labtools/NP/Run6/update_kidpar_ref.pro
;;-----------------------------------------------------------------------------------------

ps       = 1
png      = 0
charsize = 2

;; Run6 focal planes
;; kidpar_1mm_infile = !nika.off_proc_dir+"/kidpar_ref_434pixref_1mm_bestscans_v4.fits"
;; kidpar_2mm_infile = !nika.off_proc_dir+"/kidpar_ref_434pixref_2mm_bestscans_v4.fits"

;kidpar_1mm_infile = !nika.off_proc_dir+"/kidpar_ref_nopixref_1mm_bestscans_v4.fits"
;kidpar_2mm_infile = !nika.off_proc_dir+"/kidpar_ref_nopixref_2mm_bestscans_v4.fits"

;; Get new beams for Run5
run = 5
kidpar_1mm_infile = "run5_20121123s197_1mm.fits"
kidpar_2mm_infile = "run5_20121123s197_2mm.fits"
xrange_1mm = [-50,50]
yrange_1mm = [-50,60]
xrange_2mm = [-1,1]*80
yrange_2mm = [-1,1]*80

;; Get new beams for Run6
run = 6
kidpar_1mm_infile = "my_new_ref_run6_1mm.fits"
kidpar_2mm_infile = "my_new_ref_run6_2mm.fits"
xrange_1mm = [-70,70]
yrange_1mm = [-70,60]
xrange_2mm = [-1,1]*80
yrange_2mm = [-1,1]*80

;;;; Get new beams for Run Cryo
;;run = 6.5
;;kidpar_1mm_infile = !nika.off_proc_dir+"/kidpar_ref_1mm_runcryo.fits"
;;kidpar_2mm_infile = !nika.off_proc_dir+"/kidpar_ref_2mm_runcryo.fits"
;;xrange_1mm = [-70,70]
;;yrange_1mm = [-70,60]
;;xrange_2mm = [-1,1]*80
;;yrange_2mm = [-1,1]*80

;;; Keep only selected kids
;k1_ref = mrdfits(!nika.off_proc_dir+"/kidpar_ref_nopixref_1mm_bestscans_v4.fits", 1)
;k2_ref = mrdfits(!nika.off_proc_dir+"/kidpar_ref_nopixref_2mm_bestscans_v4.fits", 1)

k1 = mrdfits( kidpar_1mm_infile, 1)
k2 = mrdfits( kidpar_2mm_infile, 1)

;k1.type = k1_ref.type
;k2.type = k2_ref.type

w11 = where( k1.type eq 1, nw11)
w12 = where( k2.type eq 1, nw12)

print, nw11, nw12


phi    = dindgen(100)/99.*2*!dpi
cosphi = cos(phi)
sinphi = sin(phi)

;;beam_scale = 0.25
beam_scale = !fwhm2sigma/2 ; divide by 2 to get radius and not diameter

if ps eq 0 then wind, 1, 1, /free
;;if ps eq 1 or png eq 1 then outplot, file='run6_fpg_1mm', charsize=charsize, ps=ps, png=png
if ps eq 1 or png eq 1 then outplot, file='run'+strtrim(run,2)+'_fpg_1mm', charsize=charsize, ps=ps, png=png, thick=2, charthick=2
plot,  k1[w11].nas_x, k1[w11].nas_y, psym=1, xra=xrange_1mm, yra=yrange_1mm, /iso, /xs, /ys, $
       xtitle='Arcsec', ytitle='Arcsec', charsize=2
oplot,  k1[w11].nas_x, k1[w11].nas_y, psym=1, col=70
for i=0, nw11-1 do begin
   ikid = w11[i]

   xx1 = k1[ikid].fwhm_x*beam_scale*cosphi
   yy1 = k1[ikid].fwhm_y*beam_scale*sinphi

   x1  =  cos(k1[ikid].theta)*xx1 - sin(k1[ikid].theta)*yy1
   y1  =  sin(k1[ikid].theta)*xx1 + cos(k1[ikid].theta)*yy1

   oplot, k1[ikid].nas_x + x1, k1[ikid].nas_y + y1, col=70
endfor
legendastro, ['Run '+strtrim(long(run),2), $
              '1mm'], chars=2, /right, box=0
outplot, /close

if ps eq 0 then wind, 2, 2, /free
;;if ps eq 1 or png eq 1 then outplot, file='run6_fpg_2mm', charsize=charsize, ps=ps, png=png
if ps eq 1 or png eq 1 then outplot, file='run'+strtrim(run,2)+'_fpg_2mm', charsize=charsize, ps=ps, png=png, thick=2, charthick=2
plot,  k2[w12].nas_x, k2[w12].nas_y, psym=1, xra=xrange_2mm, yra=yrange_2mm, /iso, /xs, /ys, $
       xtitle='Arcsec', ytitle='Arcsec', charsize=2
oplot,  k2[w12].nas_x, k2[w12].nas_y, psym=1, col=250
for i=0, nw12-1 do begin
   ikid = w12[i]

   xx1 = k2[ikid].fwhm_x*beam_scale*cosphi
   yy1 = k2[ikid].fwhm_y*beam_scale*sinphi

   x1  =  cos(k2[ikid].theta)*xx1 - sin(k2[ikid].theta)*yy1
   y1  =  sin(k2[ikid].theta)*xx1 + cos(k2[ikid].theta)*yy1

   oplot, k2[ikid].nas_x + x1, k2[ikid].nas_y + y1, col=250
endfor
legendastro, ['Run '+strtrim(long(run),2), $
              '2mm'], chars=2, /right, box=0
outplot, /close

end
