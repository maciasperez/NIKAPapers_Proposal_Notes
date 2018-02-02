
whoami = 'nico' ; 'aina'

;; Outplot format
ps = 0
png= 0
nside = 256 ; 1024 ; 512

if strupcase(whoami) eq 'AINA' then begin
   input_cl_file = "/Users/andriaai/NIKA/Processing/idl-libs/Healpix/test/cl.fits"
   map_dir = "/Users/andriaai/NIKA/Processing/Labtools/AA/Dev"
endif else begin
   input_cl_file = !nika.soft_dir+"/idl-libs/Healpix/test/cl.fits"
   map_dir = "/Data/Planck/Dust"

   ;; To have B modes in the CMB simulation
   ;;readcol, '/Data/CLS_totCls.dat',  l, cl1, cl2, cl3, cl4, format='DDDDD'
   readcol, "/Users/ponthien/Projects/NIKA/Soft/Papers_Proposal_Notes/SimuKid/Pro/cmb_totCls_r0.001.dat", $
            l, cl1, cl2, cl3, cl4, format='DDDDD'
   wind, 1, 1, /free
   plot, l, cl1, /xs, /ys, /xlog, yra=[1d-6,1d5], /ylog
   oplot, l, cl2, col=70
   oplot, l, cl3, col=250
   oplot, l, abs(cl4), col=150

   ;; Need to add monopole and dipole for synfast
   cl_in = dblarr(max(l)+1,4)
   ;; prefactor included in the CAMB's output files
   fl = l*(l+1)/(2*!dpi)
   cl_in[2:*,0] = cl1/fl
   cl_in[2:*,1] = cl2/fl
   cl_in[2:*,2] = cl3/fl
   cl_in[2:*,3] = cl4/fl

   input_cl_file = 'mycl.fits'
   cl2fits, cl_in, input_cl_file
endelse

npix = nside2npix( nside)

;; Theoretical Cl in microK^2
cl_in = mrdfits( input_cl_file, 1, cl_header)
ncl_in = n_elements(cl_in)

;; Produce maps
isynfast, input_cl_file, cmb_maps, nlmax=3*nside-1, nside=nside
;for i=0, 2 do mollview, cmb_maps[*,i]
;stop

;; Tolerance on BB contamination
l_in      = dindgen(ncl_in)
fl_in     = l_in*(l_in+1)/(2*!dpi)

;; Generate spurious signals 
map_x      = dblarr(npix,3)
map_x[*,0] = cmb_maps[*,0]^2
map_x[*,1] = 2.d0*cmb_maps[*,0]*cmb_maps[*,1]
map_x[*,2] = 2.d0*cmb_maps[*,0]*cmb_maps[*,2]

ianafast, map_x, cl_x, /double, nlmax=3*nside-1, /polar, /ring
lx = dindgen(n_elements(cl_x[*,0]))
lx = lx[2:*]
clxt  = reform( cl_x[2:*,0])
clxe  = reform( cl_x[2:*,1])
clxb  = reform( cl_x[2:*,2])
clxte = reform( cl_x[2:*,3])
clxtb = reform( cl_x[2:*,4])
clxeb = reform( cl_x[2:*,5])
ncl_x = n_elements(lx)
fl_x  = lx*(lx+1)/(2.*!dpi)

;; All spectra on the same plot
yra = [1d-7, 1d9]
xra = [2, 3000]
thick=2
xra = [1, 5000]
cmb_col = !p.color
xcol = [70, 150, 250, 100, 200, 40]

if ps eq 0 then wind, 1, 1, /free, /large
outplot, file='cmb_power_spectra', ps=ps, png=png, thick=thick
plot, lx, fl_x*clxt, /xs, /ys, xra=xra, yra=yra, $
      xtitle='!12l', $
      ytitle='!12l(l+1)C!dl!n /2!7p!3 !7l!3K!u2!n', $
      /ylog, /xlog
oplot, l_in, fl_in*cl_in.temperature, col=0
oplot, lx, fl_x*clxt, col=xcol[0]
oplot, l_in, fl_in*cl_in.gradient, col=cmb_col
oplot, l_in, fl_in*cl_in.curl, col=cmb_col
oplot, l_in, fl_in*abs(cl_in.G_T), col=cmb_col
oplot, lx, fl_x*clxe, col=xcol[1]
oplot, lx, fl_x*clxb, col=xcol[2]
oplot, lx, fl_x*abs(clxte), col=xcol[3]
xyouts, 2.2, (fl_in*cl_in.temperature)[2]*1.3, 'TT'
xyouts, 2.2, (fl_in*cl_in.gradient)[2]*1.3, 'EE'
xyouts, 2.2, (fl_in*cl_in.curl)[2]*1.3, 'BB (r=0.001)'
xyouts, 2.2, (fl_in*abs(cl_in.g_t))[2]*1.3, 'TE'
legendastro, ['CMB', $
              '!7D!3TT', $
              '!7D!3EE', $
              '!7D!3BB', $
              '!7D!3TE'], line=0, $
             col=[cmb_col, xcol[0:3]]
xyouts, 10, 1d5, 'only cmb', orient=45, chars=2
outplot, /close, /verb

;;------------------------------------------------------------------------------------------
;; Look at the actual correlations between CMB and foregrounds
planck_cmb = mrdfits("/Data/Planck/Maps/COM_CMB_IQU-commander_1024_R2.02_full.fits", 1, hcmb)

;; K to microK
planck_cmb.i_stokes *= 1d6
planck_cmb.q_stokes *= 1d6
planck_cmb.u_stokes *= 1d6

;; Dust at 353GHz
planck_dust = mrdfits("/Data/Planck/Maps/COM_CompMap_DustPol-commander_1024_R2.00.fits", 1, hdust)
print, sxpar(hdust,"NU_REF")
lambda = !const.c/353.d9 * 1d6 ; microns
planck_dust.q_ml_full /= 1000.d0
planck_dust.u_ml_full /= 1000.d0
convert_millik_megajy, lambda, planck_dust.q_ml_full, q_megajy, /rj
convert_megajy_millik, lambda, q_megajy, q_dust_ref, /cmb_th
convert_millik_megajy, lambda, planck_dust.u_ml_full, u_megajy, /rj
convert_megajy_millik, lambda, u_megajy, u_dust_ref, /cmb_th
q_dust_ref *= 1000 ; mK to uK
u_dust_ref *= 1000 ; mK to uK
ud_grade, q_dust_ref, q_dust_ref_out, nside_out=nside, order_in='nested', order_out='ring'
ud_grade, u_dust_ref, u_dust_ref_out, nside_out=nside, order_in='nested', order_out='ring'
q_dust_ref = q_dust_ref_out
u_dust_ref = u_dust_ref_out

planck_t_dust = mrdfits("/Data/Planck/Maps/COM_CompMap_ThermalDust-commander_2048_R2.00.fits", 1, htdust)
print, sxpar( htdust, "NU_REF")
lambda = !const.c/545.d9 * 1d6 ; microns
mollview, planck_t_dust.i_ml_full, /nest
ud_grade, planck_t_dust.i_ml_full, i_dust_ref, nside_out=nside, order_in='nested', order_out='ring'
convert_millik_megajy, lambda, i_dust_ref, i_megajy, /rj
convert_megajy_millik, lambda, i_megajy, i_dust_ref, /cmb_th
i_dust_ref *= 1000 ; mK to uK

mollview, i_dust_ref, title='I dust (545 GHz)'
mollview, q_dust_ref, title='Q dust (353 GHz)'
mollview, u_dust_ref, title='U dust (353 GHz)'

;; Synchrotron
planck_sync = mrdfits("/Data/Planck/Maps/COM_CompMap_SynchrotronPol-commander_0256_R2.00.fits", 1, hsync)
print, sxpar( hsync, "NU_REF")
lambda = !const.c/30.d9 * 1d6 ; microns

;; uK_RJ to mK_RJ
planck_sync.q_ml_full /= 1000.d0 ; uK to mK
planck_sync.u_ml_full /= 1000.d0
convert_millik_megajy, lambda, planck_sync.q_ml_full, q_megajy, /rj
convert_megajy_millik, lambda, q_megajy, q_sync_ref, /cmb_th
convert_millik_megajy, lambda, planck_sync.u_ml_full, u_megajy, /rj
convert_megajy_millik, lambda, u_megajy, u_sync_ref, /cmb_th
q_sync_ref *= 1000 ; mK to uK
u_sync_ref *= 1000 ; mK to uK
ud_grade, q_sync_ref, q_sync_ref_out, nside_out=nside, order_in='nested', order_out='ring'
ud_grade, u_sync_ref, u_sync_ref_out, nside_out=nside, order_in='nested', order_out='ring'
q_sync_ref = q_sync_ref_out
u_sync_ref = u_sync_ref_out

planck_t_sync = mrdfits("/Data/Planck/Maps/COM_CompMap_Synchrotron-commander_0256_R2.00.fits", 1, htsync)
print, sxpar( htsync, "NU_REF")
lambda = !const.c/408.d6 * 1d6 ; microns
ud_grade, planck_t_sync.i_ml, i_sync_ref, nside_out=nside, order_in='nested', order_out='ring'
i_sync_ref /= 1000.d0 ; uK to mK
convert_millik_megajy, lambda, i_sync_ref, i_megajy, /rj
convert_megajy_millik, lambda, i_megajy, i_sync_ref, /cmb_th
i_sync_ref *= 1000 ; mK to uK

mollview, i_sync_ref, title='I sync_ref (408 MHz)'
mollview, q_sync_ref, title='Q sync_ref ( 30 GHz)'
mollview, u_sync_ref, title='U sync_ref ( 30 GHz)'

;; Cross component maps !
;; amplitude of the foreground residuals
xmap = dblarr(npix,3)

;; Mask out the galactic plane
latitude_cut = 30 ; deg
cover = dblarr(npix)
ipring = lindgen(npix)
pix2ang_ring, nside, ipring, theta, phi
latitude = 90.-theta*!radeg
w = where( abs(latitude) ge latitude_cut, nw)
cover[w] = 1.d0
fsky = float(nw)/npix

;; Cosmic Variance on BB
cosmic_var = sqrt(2./(fsky*(2*l_in+1.d0)))*cl_in.curl

;; Spurious maps
wd, /a
eta_pow = [-2, -3, -4, -5, -6]
eta_list = 10.d0^eta_pow
nu_list = [70.d0, 100.d0, 143.d0, 217.d0, 353.d0]

;eta_list = eta_list[0:1]
;nu_list = [100., 200.]
n_eta = n_elements(eta_list)
n_nu  = n_elements(nu_list)

epsilon_list = dblarr(n_eta,4)

beta_dust = 1.5
beta_sync = -3
;; ianafast, xmap, cl_x, /double, nlmax=3*nside-1, /polar, /ring
;; lx = dindgen(n_elements(cl_x[*,0]))
;; lx = lx[2:*]
;; clxt  = reform( cl_x[2:*,0])
;; clxe  = reform( cl_x[2:*,1])
;; clxb  = reform( cl_x[2:*,2])
;; clxte = reform( cl_x[2:*,3])
;; clxtb = reform( cl_x[2:*,4])
;; clxeb = reform( cl_x[2:*,5])

ps = 1
make_ct, n_eta, mycol
results_file = 'fg_res_nl.dat'
openw, lu, results_file, /get_lun
printf, lu, "#nu (GHz), eta, epsilonT, epsilonE, epsilonB, epsilonTE"
for inu=0, n_nu-1 do begin
;   outplot, file='fg_res_spectra_nu'+strtrim( long(nu_list[inu]),2), png=png, ps=ps
   plot_init = 0

   ;; Scale foregrounds to the reference frequencies
   i_dust = i_dust_ref * (nu_list[inu]/545.d0)^beta_dust
   q_dust = q_dust_ref * (nu_list[inu]/353.d0)^beta_dust
   u_dust = u_dust_ref * (nu_list[inu]/353.d0)^beta_dust

   i_sync = i_sync_ref * (nu_list[inu]/0.408)^beta_sync
   q_sync = q_sync_ref * (nu_list[inu]/30.d0)^beta_sync
   u_sync = u_sync_ref * (nu_list[inu]/30.d0)^beta_sync

   ;; Loop on residuals amplitude
   for i_eta=0, n_eta-1 do begin
      outplot, file='fg_res_spectra_nu'+strtrim( long(nu_list[inu]),2)+$
               'eta10'+strtrim(eta_pow[i_eta],2), png=png, ps=ps
      eta = eta_list[i_eta]
      
      xmap[*,0] = cmb_maps[*,0]^2 + $
                  2.*cmb_maps[*,0]*eta*(i_dust + i_sync) + $
                  eta^2*(i_dust + i_sync)^2
      xmap[*,1] = 2.d0*(cmb_maps[*,0]*cmb_maps[*,1] + $
                        eta*(q_dust+q_sync)*cmb_maps[*,0] + $
                        eta*(i_dust+i_sync)*cmb_maps[*,1])
      xmap[*,2] = 2.d0*(cmb_maps[*,0]*cmb_maps[*,2] + $
                        eta*(u_dust+u_sync)*cmb_maps[*,0] + $
                        eta*(i_dust+i_sync)*cmb_maps[*,2])

      ispice, xmap, cover, lx, clxt, clxe, clxb, clxte
      lx    = lx[2:*]
      clxt  = clxt[2:*]
      clxe  = clxe[2:*]
      clxb  = clxb[2:*]
      clxte = clxte[2:*]
      
;; Pixel window functions
      wlpix = mrdfits( !nika.soft_dir+"/idl-libs/Healpix/data/pixel_window_n"+$
                       string( nside, format='(I4.4)')+".fits", 1)
;; remove monopole and dipole
      wlpix = wlpix[2:ncl_x-1]
      
;; Derive spec
      w100 = where( l_in eq 100)
      wx100 = where( lx eq 100)
;      tol = 0.1
      
      if plot_init eq 0 then begin
;; All spectra on the same plot
         yra = [1d-7, 1d15]
         xra = [2, 3000]
         thick=2
         if ps eq 0 then wind, 1, 1, /free, /large
         plot, lx, fl_x*clxt, /xs, /ys, xra=xra, yra=yra, $
               xtitle='!12l', $
               ytitle='!12l(l+1)C!dl!n /2!7p!3 !7l!3K!u2!n', $
               /ylog, /xlog, /nodata
         xyouts, 2.2, (fl_in*cl_in.temperature)[2]*1.3, 'TT'
         xyouts, 2.2, (fl_in*cl_in.gradient)[2]*1.3, 'EE'
         xyouts, 2.2, (fl_in*cl_in.curl)[2]*1.3, 'BB (r=0.001)'
         xyouts, 2.2, (fl_in*abs(cl_in.g_t))[2]*1.3, 'TE'
;         oplot, lx, fl_x*clxt, col=xcol[0]
         oplot, l_in, fl_in*cl_in.temperature, col=0
         oplot, l_in, fl_in*cl_in.gradient, col=cmb_col
         oplot, l_in, fl_in*cl_in.curl, col=cmb_col
         oplot, l_in, fl_in*abs(cl_in.G_T), col=cmb_col
;         oplot, lx, fl_x*clxe, col=xcol[1]
;         oplot, lx, fl_x*abs(clxte), col=xcol[3]
         legendastro, ['CMB', $
                       '!7D!3TT', $
                       '!7D!3EE', $
                       '!7D!3BB', $
                       '!7D!3TE'], line=0, $
                      col=[cmb_col, xcol[0:3]]
         plot_init = 1
      endif
                                ;oplot, lx, fl_x*clxb, col=xcol[2]
      oplot, lx, fl_x*clxb, col=70 ; col=mycol[i_eta]
      w = where( lx eq 1000)
      xyouts, 1000, fl_x[w]*clxb[w]*1.1, "!7g!3 = 10!u"+strtrim( long(eta_pow),2)+"!n", col=70;, col=mycol[i_eta]

      wtol = where( l_in eq 90)
      ;; epsilon_list[i_eta] = sqrt( tol*cl_in[wtol].curl/clxb[w])
      ;; epsilon_list[i_eta] = sqrt( cosmic_var[wtol]/clxb[w])

      cosmic_var = sqrt(2./(fsky*(2*l_in+1.d0)))*cl_in.temperature
      epsilon_list[i_eta,0] = sqrt( cosmic_var[wtol]/clxt[w])

      cosmic_var = sqrt(2./(fsky*(2*l_in+1.d0)))*cl_in.gradient
      epsilon_list[i_eta,1] = sqrt( cosmic_var[wtol]/clxe[w])

      cosmic_var = sqrt(2./(fsky*(2*l_in+1.d0)))*cl_in.curl
      epsilon_list[i_eta,2] = sqrt( cosmic_var[wtol]/clxb[w])

      cosmic_var = sqrt(1./(fsky*(2*l_in+1.d0)))*sqrt(cl_in.g_t^2 + cl_in.temperature*cl_in.gradient)
      epsilon_list[i_eta,3] = sqrt( cosmic_var[wtol]/abs(clxte[w]))

      ;; Output ascii file
      printf, lu, strtrim(nu_list[inu],2)+", 10^"+$
              strtrim(eta_pow[i_eta],2)+", "+$
              strtrim(epsilon_list[i_eta,0],2)+", "+$
              strtrim(epsilon_list[i_eta,1],2)+", "+$
              strtrim(epsilon_list[i_eta,2],2)+", "+$
              strtrim(epsilon_list[i_eta,3],2)


;;;for itol=0, ntol-1 do oplot, lx, epsilon[itol]^2*fl_x*clxb, col=200
;;legendastro, '!7e!3!dres!n = '+strtrim(eta,2)]
;;for itol=0, ntol-1 do leg_txt = [leg_txt, "!7e!3!dspec!n: "+
;;                             strtrim(eta,2)+", "+$
;;                             strtrim(epsilon[itol],2)
;;

      outplot, /close, /verb
   endfor
;   legendastro, [strtrim(long(nu_list[inu]),2)+' GHz', $
;                 "!7g!3 = 10!u"+strtrim( long(eta_pow),2)+"!n, !7e!3 = "+strtrim(epsilon_list,2)], /right
;   outplot, /close, /verb
endfor
close, lu
free_lun, lu
spawn, "cat "+results_file

readcol, results_file, nu, eta_test, epst, epse, epsb, epste, $
         format='A,A,A,A,A,A', delim=',', comment='#'
all_eps = dblarr(4,n_elements(nu))
all_eps[0,*] = double(epst)
all_eps[1,*] = double(epse)
all_eps[2,*] = double(epsb)
all_eps[3,*] = double(epste)
spec_title = ['T', 'E', 'B', 'TE']
make_ct, n_nu, ct
if ps eq 0 then wind, 1, 1, /free, /large
my_multiplot, 2, 2, pp, pp1, /rev
outplot, file='fg_res_nl_plot', ps=ps, png=png
xra = minmax(eta_list)*[0.1,10]
yra = minmax(all_eps)
for ispec=0, 3 do begin
   for i_nu=0, n_nu-1 do begin
      w = where( double(nu) eq nu_list[i_nu], nw)
      if i_nu eq 0 then begin
         plot, eta_list, all_eps[ispec,w], /xlog, /ylog, xra=xra, yra=yra, $
               xtitle='!7g!3', ytitle='!7e!3 ('+strtrim(spec_title[ispec],2)+')', $
               position=pp1[ispec,*], /noerase
         legendastro, strtrim( long(nu_list),2)+" GHz", col=ct, line=0, /right
      endif
      oplot, eta_list, all_eps[ispec,w], psym=-8, col=ct[i_nu]
   endfor
endfor
outplot, /close, /verb

end
