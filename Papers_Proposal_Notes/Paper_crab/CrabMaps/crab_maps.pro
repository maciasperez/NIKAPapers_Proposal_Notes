;; XPOL

dir_xpol = '/Users/ritacco/Desktop/maps/'

map_ja_i = mrdfits(dir_xpol+'crab_I_coadd28_reproj_nbl_ccentonly.fits',0,head_ja)
map_ja_q = mrdfits(dir_xpol+'crab_Q_coadd28_reproj_nbl_ccentonly.fits',0,head_ja)
map_ja_u = mrdfits(dir_xpol+'crab_U_coadd28_reproj_nbl_ccentonly.fits',0,head_ja)
extast,head_ja,astr
reso_ja = astr.cdelt[1]*3600.0


;; SCUPOL
dir_scuba = '/Users/ritacco/SCUBAPOL/'
map_scuba_i = mrdfits(dir_scuba+'scupollegacy_crab_cube.fits', 0, head_scuba)
extast, head_scuba, astr
reso_scuba = astr.cdelt[1]*3600.0


dirmaps = "/Users/macias/NIKAdocuments/CrabPolar/"

; Helene maps

hi = mrdfits(dirmaps+"crab_all_ptot_nika1mm_0.fits", 0,header)

wbad = where(finite(hi[*,*,0]) eq 0 or finite(hi[*,*,1]) eq 0,nwbad)
 
ih = reform(hi[*,*,0])
ih[wbad]= 0.0
ihv = reform(hi[*,*,1])
ihv[wbad]=0

extast, header, hast
; Alessia's maps


m1mm = mrdfits(dirmaps+"MAPS_1_Crab_v3.fits",1,h)
extast, h, ast
ast.ctype = hast.ctype
3ev1mm = mrdfits(dirmaps+"MAPS_1_Crab_v3.fits",2,h)
h1mm = mrdfits(dirmaps+"MAPS_1_Crab_v3.fits",3,h)


; work in astrometry



nx = ast.naxis[0]
ny = ast.naxis[1]
xarr = lindgen(nx)# replicate(1,ny)
yarr = replicate(1,nx) # lindgen(ny) 
XY2AD, xarr, yarr, ast, ra, dec
AD2XY,ra,dec,hast,hxarr,hyarr

hnx = hast.naxis[0]
hny = hast.naxis[1]

mnq = dblarr(hnx,hny)
hnq = dblarr(hnx,hny)
vnq = dblarr(hnx,hny)
hnu = dblarr(hnx,hny)
mnu = dblarr(hnx,hny)
vnu = dblarr(hnx,hny)

for ix=0,nx-1 do begin
   for iy=0,ny-1 do begin
      if v1mm[xarr[ix,iy],yarr[ix,iy],1] gt 0.0 then begin 
         mnq[hxarr[ix,iy],hyarr[ix,iy]] += m1mm[xarr[ix,iy],yarr[ix,iy],1]/v1mm[xarr[ix,iy],yarr[ix,iy],1]/v1mm[xarr[ix,iy],yarr[ix,iy],1]
         hnq[hxarr[ix,iy],hyarr[ix,iy]] += (1/v1mm[xarr[ix,iy],yarr[ix,iy],1]/v1mm[xarr[ix,iy],yarr[ix,iy],1])
      endif

      if  v1mm[xarr[ix,iy],yarr[ix,iy],2] gt 0.0 then begin 
         mnu[hxarr[ix,iy],hyarr[ix,iy]] += m1mm[xarr[ix,iy],yarr[ix,iy],2]/v1mm[xarr[ix,iy],yarr[ix,iy],2]/v1mm[xarr[ix,iy],yarr[ix,iy],2]
         hnu[hxarr[ix,iy],hyarr[ix,iy]] += (1/v1mm[xarr[ix,iy],yarr[ix,iy],2]/v1mm[xarr[ix,iy],yarr[ix,iy],2])  
      endif 
   endfor
endfor
w = where(hnq gt 0,nw)
mnq[w] /= hnq[w]
vnq[w] = sqrt(1.0/hnq[w])
w = where(hnu gt 0,nw)
mnu[w] /= hnu[w]
vnu[w] = sqrt(1.0/hnu[w])

dispim_bar, mnq, /nocon, /asp,cr=[-0.1,0.1]
dispim_bar, mnu, /nocon, /asp,cr=[-0.1,0.1]


map = dblarr(hnx,hny,3)
stdmap = dblarr(hnx,hny,3)
map[*,*,0] = ih
stdmap[*,*,0] = ihv

map[*,*,1] = mnq
map[*,*,2] = mnu
stdmap[*,*,1] = vnq
stdmap[*,*,2] = vnu




mwrfits,0, "ale2h_Crab_1mm.fits", primaryheader, /create
mwrfits,map, "ale2h_Crab_1mm.fits", header
mwrfits,stdmap, "ale2h_Crab_1mm.fits", header



; 2 mm
m2mm = mrdfits(dirmaps+"MAPS_2_Crab_v3.fits",1,h)
ast.ctype = hast.ctype
v2mm = mrdfits(dirmaps+"MAPS_2_Crab_v3.fits",2,h)
h2mm = mrdfits(dirmaps+"MAPS_2_Crab_v3.fits",3,h)

mni = dblarr(hnx,hny)
hni = dblarr(hnx,hny)
vni = dblarr(hnx,hny)
mnq = dblarr(hnx,hny)
hnq = dblarr(hnx,hny)
vnq = dblarr(hnx,hny)
hnu = dblarr(hnx,hny)
mnu = dblarr(hnx,hny)
vnu = dblarr(hnx,hny)

for ix=0,nx-1 do begin
   for iy=0,ny-1 do begin
      if v2mm[xarr[ix,iy],yarr[ix,iy],0] gt 0.0 then begin 
         mni[hxarr[ix,iy],hyarr[ix,iy]] += m2mm[xarr[ix,iy],yarr[ix,iy],0]/v2mm[xarr[ix,iy],yarr[ix,iy],0]/v2mm[xarr[ix,iy],yarr[ix,iy],0]
         hni[hxarr[ix,iy],hyarr[ix,iy]] += (1/v2mm[xarr[ix,iy],yarr[ix,iy],0]/v2mm[xarr[ix,iy],yarr[ix,iy],0])
      endif
      if v2mm[xarr[ix,iy],yarr[ix,iy],1] gt 0.0 then begin 
         mnq[hxarr[ix,iy],hyarr[ix,iy]] += m2mm[xarr[ix,iy],yarr[ix,iy],1]/v2mm[xarr[ix,iy],yarr[ix,iy],1]/v2mm[xarr[ix,iy],yarr[ix,iy],1]
         hnq[hxarr[ix,iy],hyarr[ix,iy]] += (1/v2mm[xarr[ix,iy],yarr[ix,iy],1]/v2mm[xarr[ix,iy],yarr[ix,iy],1])
      endif

      if  v2mm[xarr[ix,iy],yarr[ix,iy],2] gt 0.0 then begin 
         mnu[hxarr[ix,iy],hyarr[ix,iy]] += m2mm[xarr[ix,iy],yarr[ix,iy],2]/v2mm[xarr[ix,iy],yarr[ix,iy],2]/v2mm[xarr[ix,iy],yarr[ix,iy],2]
         hnu[hxarr[ix,iy],hyarr[ix,iy]] += (1/v2mm[xarr[ix,iy],yarr[ix,iy],2]/v2mm[xarr[ix,iy],yarr[ix,iy],2])  
      endif 
   endfor
endfor

w = where(hni gt 0,nw)
mni[w] /= hni[w]
vni[w] = sqrt(1.0/hni[w])
w = where(hnq gt 0,nw)
mnq[w] /= hnq[w]
vnq[w] = sqrt(1.0/hnq[w])
w = where(hnu gt 0,nw)
mnu[w] /= hnu[w]
vnu[w] = sqrt(1.0/hnu[w])

dispim_bar, mnq, /nocon, /asp,cr=[-0.1,0.1]
dispim_bar, mnu, /nocon, /asp,cr=[-0.1,0.1]


map = dblarr(hnx,hny,3)
stdmap = dblarr(hnx,hny,3)

map[*,*,0] = mni
map[*,*,1] = mnq
map[*,*,2] = mnu
stdmap[*,*,0] = vni
stdmap[*,*,1] = vnq
stdmap[*,*,2] = vnu


mwrfits,0, "ale2h_Crab_2mm.fits", primaryheader, /create
mwrfits,map, "ale2h_Crab_2mm.fits", header
mwrfits,stdmap, "ale2h_Crab_2mm.fits", header

