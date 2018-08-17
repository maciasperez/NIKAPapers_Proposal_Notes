
readcol, '$SK_DIR/Cl/cmb_totCls_r0.001.dat', l, clt, cle, clb, clte, $
         format='D,D,D,D', comment='#'

nside = 256

;; Need to add monopole and dipole for synfast
cl_in = dblarr(max(l)+1,4)
;; prefactor included in the CAMB's output files
fl = l*(l+1)/(2*!dpi)
cl_in[2:*,0] = clt/fl
cl_in[2:*,1] = cle/fl
cl_in[2:*,2] = clb/fl
cl_in[2:*,3] = clte/fl

input_cl_file = 'mycl.fits'
cl2fits, cl_in, input_cl_file

npix = nside2npix( nside)
cl_in = mrdfits( input_cl_file, 1, cl_header)
ncl_in = n_elements(cl_in)
isynfast, input_cl_file, cmb_maps, nlmax=3*nside-1, nside=nside

;; no mask (OK)
cover = dblarr(npix)+1.d0

;; Mask out the galactic plane (OK too)
latitude_cut = 30 ; deg
cover        = dblarr(npix)
ipring       = lindgen(npix)
pix2ang_ring, nside, ipring, theta, phi
latitude = 90.-theta*!radeg
w = where( abs(latitude) ge latitude_cut, nw)
cover[w] = 1.d0

ispice, cmb_maps, cover, l_out, clt_out, cle_out, clb_out, clte_out
l_out    = l_out[   2:*]
clt_out  = clt_out[ 2:*]
cle_out  = cle_out[ 2:*]
clb_out  = clb_out[ 2:*]
clte_out = clte_out[2:*]

;; Correct for Healpix's pixel window
;; function                                                                                             
wl = mrdfits( !healpix.path.data+"/pixel_window_n"+string(nside,form='(I4.4)')+".fits",1)
clt_out  = clt_out[ 2:*]/wl.temperature^2
cle_out  = cle_out[ 2:*]/wl.polarization^2
clb_out  = clb_out[ 2:*]/wl.polarization^2
clte_out = clte_out[2:*]/wl.temperature/wl.polarization

fl = l_out*(l_out+1)/(2.*!dpi)

xra = [1,3*nside]
yra = [1d-8,1d4]
col_ee = 70
col_bb = 250
col_te = 150
wind, 1, 1, /free, /large
plot, xra, yra, xtitle='!12l!3', ytitle='!12l(l+1)C!dl!n/2!7p!3', $
      /nodata, /xlog, /ylog, /xs, /ys
oplot, l, clt
oplot, l, cle
oplot, l, clb
oplot, l, abs(clte)
oplot, l_out, fl*clt_out
oplot, l_out, fl*cle_out, col=col_ee
oplot, l_out, fl*clb_out, col=col_bb
oplot, l_out, fl*abs(clte_out), col=col_te
legendastro, ['T', 'E', 'B', 'TE'], col=[!p.color, col_ee, col_bb, col_te], line=0

;; Check numerical precision
nmc = 30
for imc=0, nmc-1 do begin
   isynfast, input_cl_file, cmb_maps, nlmax=3*nside-1, nside=nside
   ispice, cmb_maps, cover, l_out, clt_out, cle_out, clb_out, clte_out
   l_out    = l_out[   2:*]
   clt_out  = clt_out[ 2:*]
   cle_out  = cle_out[ 2:*]
   clb_out  = clb_out[ 2:*]
   clte_out = clte_out[2:*]
   
;; Correct for Healpix's pixel window function
   wl = mrdfits( !healpix.path.data+"/pixel_window_n"+string(nside,form='(I4.4)')+".fits",1)
   clt_out  = clt_out[ 2:*]/wl.temperature^2
   cle_out  = cle_out[ 2:*]/wl.polarization^2
   clb_out  = clb_out[ 2:*]/wl.polarization^2
   clte_out = clte_out[2:*]/wl.temperature/wl.polarization

   if imc eq 0 then begin
      clt_out_res  = dblarr(nmc,n_elements(clt_out))
      cle_out_res  = clt_out_res
      clb_out_res  = clt_out_res
      clte_out_res = clt_out_res
   endif
   clt_out_res[imc,*]  = clt_out
   cle_out_res[imc,*]  = cle_out
   clb_out_res[imc,*]  = clb_out
   clte_out_res[imc,*] = clte_out
endfor
;; Plot average results
mc_reduce, clt_out_res, clt_out_avg, sigma_clt_out_avg
mc_reduce, cle_out_res, cle_out_avg, sigma_cle_out_avg
mc_reduce, clb_out_res, clb_out_avg, sigma_clb_out_avg
mc_reduce, clte_out_res, clte_out_avg, sigma_clte_out_avg
;; Final plot
fl = l_out*(l_out+1)/(2.*!dpi)
xra = [1,3*nside]
yra = [1d-8,1d4]
col_ee = 70
col_bb = 250
col_te = 150
wind, 1, 1, /free, /large
plot, xra, yra, xtitle='!12l!3', ytitle='!12l(l+1)C!dl!n/2!7p!3', $
      /nodata, /xlog, /ylog, /xs, /ys
xyouts, 10, 1d-2, 'Pure Synfast simu B=0', orient=30
oplot, l, clt
oplot, l, cle
oplot, l, clb
oplot, l, abs(clte)
oplot, l_out, fl*clt_out_avg
oplot, l_out, fl*cle_out_avg, col=col_ee
oplot, l_out, fl*clb_out_avg, col=col_bb
oplot, l_out, fl*abs(clte_out_avg), col=col_te
legendastro, ['T', 'E', 'B', 'TE'], col=[!p.color, col_ee, col_bb, col_te], line=0
legendastro, 'Nmc = '+strtrim(nmc,2), /right

end



