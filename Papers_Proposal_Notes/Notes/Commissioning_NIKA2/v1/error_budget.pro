
f = 1d9
delta = 45.d0 ; elevation
c0 = 1d9
c1 = 1d3
T0 = 270.d0

source = '0923+392'
run = 32
readcol, !nika.soft_dir+"/Labtools/AA/Polar/Quasars/ScanLists/Run"+$
         strtrim(run,2)+'/'+strupcase(source)+'_used_scan_list.txt', scan_list,$
         format='A', comment='#', /silent
scan_list = scan_list[where(scan_list ne '20180614s298' and $
                                  scan_list ne '20180614s304' and $
                                  scan_list ne '20180614s305' and $
                                  scan_list ne '20180615s296' and $
                                  scan_list ne '20180613s201' and $
                                  scan_list ne '20180614s292')]
scan_list = scan_list[ where(scan_list ne '20180616s336' and $
                                   scan_list ne '20180614s245' and $
                                   scan_list ne '20180614s340')]
scan = scan_list[0]

nk_default_param, param
nk_default_info, info
param.scan = scan
nk_getdata, param, info, data, kidpar, /polar
stop








nk_get_kidpar_ref, s, d, i, scan='20180904s1', kidpar=kidpar
w1 = where( kidpar.array eq 3 and kidpar.type eq 1 and kidpar.c0_skydip ne 0.d0, nw1)

c0 = median( kidpar[w1].c0_skydip)
c1 = median( kidpar[w1].c1_skydip)

;; "stddev"... sic !
sigma_c0 = stddev( kidpar[w1].c0_skydip)
sigma_c1 = stddev( kidpar[w1].c1_skydip)

wind, 1, 1, /free, /large
my_multiplot, 2, 2, pp, pp1, /rev
plot, kidpar[w1].c0_skydip, /xs, position=pp1[0,*]
oplot, indgen(nw1), dblarr(nw1)+c0
oplot, indgen(nw1), dblarr(nw1)+c0-sigma_c0, line=2
oplot, indgen(nw1), dblarr(nw1)+c0+sigma_c0, line=2
legendastro, 'C0'

plot, kidpar[w1].c1_skydip, /xs, position=pp1[1,*], /noerase
oplot, indgen(nw1), dblarr(nw1)+c1
oplot, indgen(nw1), dblarr(nw1)+c1-sigma_c1, line=2
oplot, indgen(nw1), dblarr(nw1)+c1+sigma_c1, line=2
legendastro, 'C1'

nmc = 100
alpha = randomn( seed, nmc)*sigma_c0
beta  = randomn( seed, nmc)*sigma_c1
;;tau = sin(delta*!dtor) * alog( 1.d0-(f-c0*(1+alpha))/(T0*c1*(1.+beta)))




np_histo, tau, position=pp1[2,*], /noerase, /fill

end
