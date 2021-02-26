;
fcut = 1./3                      ; Hz
fov = 6.5*60. ; arcsec
speednom = fov*fcut ; 130
speed = 150.
fsampling = 23.*[1, 2]           ; Hz
sample = speed/fsampling  ;8.7,4.3
print,exp(0.5/sin(45/!radeg)) ; 2
print,20.E-4/60/sqrt(0.7)     ; 10 mJys1/2 * 2(opa) * 1e-4 -->
                                ; 4E-5 y 1sigma in hr1/2 /beam

print, 20*12./(!pi*6.5^2/4)      ; 7 FOV reduce the sensitivity to
print,20.E-4/60/sqrt(0.7)*sqrt(7) ; 1.1E-4 per beam in one hour
print,20.E-4/60/sqrt(0.7)*sqrt(7)/sqrt(2) ; 7.5E-5 per beam in 2 hours
print,20.E-4/60/sqrt(0.7)*sqrt(7)/sqrt(5) ; 5E-5 per beam in 5 hours

print, 20E-3/12/sqrt(3600)*sqrt(0.7)
print, 20E-3/12/sqrt(3600)*sqrt(0.7)*sqrt(7/2.)

; Number of hours
frac = 0.25
; 3 months of NIKA2
print, 24.*365*frac* 3/12.
; 2200hrs --> 550 hrs

; Profile
beta = 0.86
print, (1-3*beta)/2  ; -0.79
Thetacore = 60.
nstep = 60
theta = 5*findgen(nstep)
y0 = 3.5E-4  ; no beam deconvolved: Reese: 1900.E-6/2.7/2.
yprof = y0*(1+ (theta/thetacore)^2)^((1-beta*3)/2.)
error = 4.3E-5                    ; for one beam
beam = 18.  ; "
;annulus integration
dthe = 30.                      ; theta difference
err_theta = dthe*findgen(nstep/6)
area = 2*!pi*err_theta*dthe
err_ann = error/sqrt(area/ (!pi/4 * beam^2))
err_ann[0] = !values.F_nan
loadct, 15
fxd_ps, /color
plot, theta, yprof/1E-4, /ylog,  ytitle = 'y [1E-4]', yra = [1e-2, 5], $
      xtitle = 'Radius [arcsec]'
oplot, err_theta, err_ann/1E-4, col = 100, psym = 4, thick = 2
oplot, theta_FK, y_FK/1E-4, col = 200
oplot, [8, 8.]*60./2, [1E-4, 10],  line = 2
fxd_psout, save = 'Arsenal_SZprofile.pdf', /rotate


; Statistics of clusters y500
cldb = READ_CSV( '/home/desertf/Desktop/ProposalIRAM_nika2/Arsenal/SZclusterDBcrop.csv', count=ncl,  header = nam, table_header = allnam) ;,  n_table_header = 1)
u = where( strmatch(nam, '*y5r500*'), nu) & print, nu
print, nam[u]                   ;psz2_y5r500 psz2_y5r500_err
y5 = cldb.(u[0])
y5e = cldb.(u[1])
u = where( strmatch(nam, '*psz2_snr*'), nu) & print, nu
snr = cldb.(u[0])
u = where( strmatch(nam, '*psz2_ra*'), nu) & print, nu
radeg = cldb.(u[0])
u = where( strmatch(nam, '*psz2_dec*'), nu) & print, nu
decdeg = cldb.(u[0])
;v = where( y5 gt 2.*y5e, nv) & print, nv  ; 309 above 3sigma, 533 above 2sigma
v = where( snr ge 3 and decdeg gt (-15.), nv) & print, nv ; 400
; Try to select the redshift range 0.15 to 0.3
histo_make, y5[v] < 13., /plot, /stat, $
            n_bins = 14, minval = 0.5, maxval = 13.5, xarr,  yarr

fxd_ps
plot, xarr, yarr, ysty = 2, ytitle = 'Number of clusters', xtitle = 'Y5R500 [1E-3 arcmin^2]', psym = 10
fxd_psout, save = 'Arsenal_SZclusters.pdf', /rotate


;; Salut Xavier,

;; Tu as utilisé un beta-model pour le profil de densité (du coup avec un modèle isotherme pour avoir la pression) ou pour le profil de Compton directement ?

;; Ci-dessous, le profil que je trouve pour un profil de pression Arnaud et al. 2010, avec les caractéristiques de A2163 dans Planck (Y_500 = 0.0131 arcmin^2, soit M_500 = 18.74 e14 Msun). Les angles sont en secondes d’arc, et le profil n’est convolué ni par le beam, ni par la fonction de transfert. Pour la conversion, c’est environ -11.9 Jy/beam/y.

;; N’hésite pas si tu as besoin de plus !

;; Bonne soirée,
;; Florian

theta_FK = [0.,   3.,   6.,   9.,  12.,  15.,  18.,  21.,  24.,  27.,  30., 33.,  36.,  39.,  42.,  45.,  48.,  51.,  54.,  57.,  60.,  63., 66.,  69.,  72.,  75.,  78.,  81.,  84.,  87.,  90.,  93.,  96., 99., 102., 105., 108., 111., 114., 117., 120., 123., 126., 129., 132., 135., 138., 141., 144., 147., 150., 153., 156., 159., 162., 165., 168., 171., 174., 177., 180., 183., 186., 189., 192., 195., 198., 201., 204., 207., 210., 213., 216., 219., 222., 225., 228., 231., 234., 237., 240., 243., 246., 249., 252., 255., 258., 261., 264., 267., 270., 273., 276., 279., 282., 285., 288., 291., 294., 297., 300.]

y_FK = [4.71701399e-04, 4.54216991e-04, 4.42342906e-04, 4.31607979e-04, 4.21487750e-04, 4.11789184e-04, 4.02420174e-04, 3.93330148e-04, 3.84488576e-04, 3.75875619e-04, 3.67477540e-04, 3.59284257e-04, 3.51287956e-04, 3.43482273e-04, 3.35861788e-04, 3.28421707e-04, 3.21157657e-04, 3.14065557e-04, 3.07141525e-04, 3.00381825e-04, 2.93782825e-04, 2.87340977e-04, 2.81052796e-04, 2.74914853e-04, 2.68923772e-04, 2.63076219e-04, 2.57368909e-04, 2.51798601e-04, 2.46362100e-04, 2.41056259e-04, 2.35877979e-04, 2.30824210e-04, 2.25891953e-04, 2.21078263e-04, 2.16380246e-04, 2.11795063e-04, 2.07319931e-04, 2.02952120e-04, 1.98688959e-04, 1.94527829e-04, 1.90466172e-04, 1.86501483e-04, 1.82631314e-04, 1.78853272e-04, 1.75165021e-04, 1.71564278e-04, 1.68048816e-04, 1.64616460e-04, 1.61265090e-04, 1.57992636e-04, 1.54797081e-04, 1.51676458e-04, 1.48628850e-04, 1.45652388e-04, 1.42745252e-04, 1.39905669e-04, 1.37131909e-04, 1.34422292e-04, 1.31775177e-04, 1.29188970e-04, 1.26662116e-04, 1.24193103e-04, 1.21780458e-04, 1.19422749e-04, 1.17118579e-04, 1.14866591e-04, 1.12665463e-04, 1.10513910e-04, 1.08410678e-04, 1.06354551e-04, 1.04344342e-04, 1.02378898e-04, 1.00457096e-04, 9.85778447e-05, 9.67400805e-05, 9.49427690e-05, 9.31849039e-05, 9.14655056e-05, 8.97836210e-05, 8.81383227e-05, 8.65287080e-05, 8.49538987e-05, 8.34130399e-05, 8.19052998e-05, 8.04298687e-05, 7.89859586e-05, 7.75728028e-05, 7.61896546e-05, 7.48357876e-05, 7.35104946e-05, 7.22130870e-05, 7.09428947e-05, 6.96992653e-05, 6.84815637e-05, 6.72891714e-05, 6.61214864e-05, 6.49779225e-05, 6.38579087e-05, 6.27608893e-05, 6.16863228e-05, 6.06336822e-05]
