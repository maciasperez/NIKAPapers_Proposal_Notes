pro nika2_sz_time_estim, y_at_r500, t500, path_t_obs=path_t_obs, small_map=small_map, sensitivity=sensitivity, snr=snr, valid_kids=0.6

;y_at_r500 in mJy
;t500 in arcmin

if not keyword_set(sensitivity) then sensitivity=20.;mJy
if not keyword_set(SNR) then snr=4 ;snr at t500
if not keyword_set(valid_kids) then valid_kids=0.6 ;snr at t500

restore,  path_t_obs+'/int_time_radecscan_37.7min_13x8.sav' ; maps of obs. time on the area (considering pixels of 2arcsec)

if keyword_set(small_map) then restore,  path_t_obs+'/int_time_radecscan_25.7min_11x7.sav'; maps of obs. time on the area (considering pixels of 2arcsec)
if keyword_set(small_map) then t_one_scan=27. else t_one_scan=40.

  fov=10.
  reso=2./60
  nmap = 2*long(10/reso/2)+1
  dd = dblarr(nmap,nmap)
  for i=0l,nmap-1 do begin      
     for j=0l,nmap-1 do begin
        rad = reso*((double(i-(nmap-1)/2.0))^2 + (double(j-(nmap-1)/2.0))^2)^0.5 
        dd[i,j] = rad
     endfor
  endfor

rs=(findgen(15)*20.+10)/60.
rs=[0,rs]
tpp_r=fltarr(n_elements(rs))
n_ind=fltarr(n_elements(rs))
  for i=0, n_elements(rs)-2 do begin
     tpp_r[i]=mean(tpp[where(DD le rs[i+1] and DD gt rs[i])])
  endfor

  at_r500=where(abs(rs-t500) eq min(abs(rs-t500)))

   pix_per_beam=!dpi*18.*18./(!dpi*2.*2) 

   itr=(sensitivity*sensitivity*snr*snr/((y_at_r500)^2.))/(valid_kids*pix_per_beam*tpp_r[at_r500[0]]) ;number of iteration of the simulated otf required
   hours=itr*t_one_scan/60.

print, strtrim(hours,2)+' hours'

end
