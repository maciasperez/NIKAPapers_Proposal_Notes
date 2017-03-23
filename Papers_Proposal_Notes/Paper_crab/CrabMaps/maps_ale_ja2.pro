;; XPOL

dir_xpol = '/Users/ritacco/Desktop/maps/'

map_ja_i = mrdfits(dir_xpol+'crab_I_coadd28_reproj_nbl_ccentonly.fits',0,head_ja)
map_ja_q = mrdfits(dir_xpol+'crab_Q_coadd28_reproj_nbl_ccentonly.fits',0,head_ja)
map_ja_u = mrdfits(dir_xpol+'crab_U_coadd28_reproj_nbl_ccentonly.fits',0,head_ja)
extast,head_ja,astr_ja
astr_ja.ctype[0] = "RA---TAN"
astr_ja.ctype[1] = "DEC--TAN"

reso_ja = astr_ja.cdelt[1]*3600.0

; work in astrometry

; 2 mm
dir_2mm = '/Volumes/STIPU/NIKA/Plots/CRAB_RADEC_LkgCorr/'
m2mm = mrdfits(dir_2mm+"MAPS_2_Crab_v3.fits",1,h)
head = h
sxaddpar, head, "naxis", 2
sxdelpar, head, "naxis3"
sxaddpar, head, "ctype1", "RA---TAN"
sxaddpar, head, "ctype2", "DEC--TAN"

head = head_ni
sxaddpar, head, "naxis", 2
sxdelpar, head, "naxis3"
sxaddpar, head, "ctype1", "RA---TAN"
sxaddpar, head, "ctype2", "DEC--TAN"
map = reform(map_2mm[*,*,0])
delta_x = -85
delta_y = -5
map = shift( map, delta_x, delta_y)

extast,h,ast
ast.ctype[0] = "RA---TAN"
ast.ctype[1] = "DEC--TAN"
ast.longpole = 0.d0
v2mm = mrdfits(dir_2mm+"MAPS_2_Crab_v3.fits",2,h)
h2mm = mrdfits(dir_2mm+"MAPS_2_Crab_v3.fits",3,h)


nx = ast.naxis[0]
ny = ast.naxis[1]

xarr = lindgen(nx)# replicate(1,ny)
yarr = replicate(1,nx) # lindgen(ny) 
XY2AD, xarr, yarr, ast, ra, dec

AD2XY,ra,dec,astr_ja,hxarr,hyarr

hnx = astr_ja.naxis[0]
hny = astr_ja.naxis[1]

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
    if hxarr[ix,iy] ge 0 and hxarr[ix,iy] lt hnx and  hyarr[ix,iy] ge 0 and hyarr[ix,iy] lt hny then begin 
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

dispim_bar, mnq, /nocon, /asp,cr=[-0.1,0.2]
dispim_bar, mnu, /nocon, /asp,cr=[-0.1,0.2]


map = dblarr(hnx,hny,3)
stdmap = dblarr(hnx,hny,3)

map[*,*,0] = mni
map[*,*,1] = mnq
map[*,*,2] = mnu
stdmap[*,*,0] = vni
stdmap[*,*,1] = vnq
stdmap[*,*,2] = vnu


mwrfits,0, "ale2h_reso_xpol_Crab_2mm.fits", primaryheader, /create
mwrfits,map, "ale2h_reso_xpol_Crab_2mm.fits", head_ja
mwrfits,stdmap, "ale2h_reso_xpol_Crab_2mm.fits", head_ja


;; SCUPOL
dir_scuba = '/Users/ritacco/SCUBAPOL/'
m850 = mrdfits(dir_scuba+'scupollegacy_crab_cube.fits', 0, h)
m850 = m850*450.d0
extast,h,ast
reso_scuba = ast.cd[3]*3600.0

ast.ctype[0] = "RA---TAN"
ast.ctype[1] = "DEC--TAN"
ast.longpole = 0.d0
v850 = mrdfits(dir_scuba+'scupollegacy_crab_cube.fits',1,h)
v850 = v850*450.d0

nx = ast.naxis[0]
ny = ast.naxis[1]

xarr = lindgen(nx)# replicate(1,ny)
yarr = replicate(1,nx) # lindgen(ny) 
XY2AD, xarr, yarr, ast, ra, dec

AD2XY,ra,dec,astr_ja,hxarr,hyarr

hnx = astr_ja.naxis[0]
hny = astr_ja.naxis[1]

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
    if hxarr[ix,iy] ge 0 and hxarr[ix,iy] lt hnx and  hyarr[ix,iy] ge 0 and hyarr[ix,iy] lt hny then begin 
      if v850[xarr[ix,iy],yarr[ix,iy],0] gt 0.0 then begin 
         mni[hxarr[ix,iy],hyarr[ix,iy]] += m850[xarr[ix,iy],yarr[ix,iy],0]/v850[xarr[ix,iy],yarr[ix,iy],0]/v850[xarr[ix,iy],yarr[ix,iy],0]
         hni[hxarr[ix,iy],hyarr[ix,iy]] += (1/v850[xarr[ix,iy],yarr[ix,iy],0]/v850[xarr[ix,iy],yarr[ix,iy],0])
      endif
      if v850[xarr[ix,iy],yarr[ix,iy],1] gt 0.0 then begin 
         mnq[hxarr[ix,iy],hyarr[ix,iy]] += m850[xarr[ix,iy],yarr[ix,iy],1]/v850[xarr[ix,iy],yarr[ix,iy],1]/v850[xarr[ix,iy],yarr[ix,iy],1]
         hnq[hxarr[ix,iy],hyarr[ix,iy]] += (1/v850[xarr[ix,iy],yarr[ix,iy],1]/v850[xarr[ix,iy],yarr[ix,iy],1])
      endif

      if  v850[xarr[ix,iy],yarr[ix,iy],2] gt 0.0 then begin 
         mnu[hxarr[ix,iy],hyarr[ix,iy]] += m850[xarr[ix,iy],yarr[ix,iy],2]/v850[xarr[ix,iy],yarr[ix,iy],2]/v850[xarr[ix,iy],yarr[ix,iy],2]
         hnu[hxarr[ix,iy],hyarr[ix,iy]] += (1/v850[xarr[ix,iy],yarr[ix,iy],2]/v850[xarr[ix,iy],yarr[ix,iy],2])  
      endif 
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

dispim_bar, mnq, /nocon, /asp,cr=[-0.1,0.2]
dispim_bar, mnu, /nocon, /asp,cr=[-0.1,0.2]


map = dblarr(hnx,hny,3)
stdmap = dblarr(hnx,hny,3)

map[*,*,0] = mni
map[*,*,1] = mnq
map[*,*,2] = mnu
stdmap[*,*,0] = vni
stdmap[*,*,1] = vnq
stdmap[*,*,2] = vnu


mwrfits,0, "ale2h_reso_xpol_Crab_850.fits", primaryheader, /create
mwrfits,map, "ale2h_reso_xpol_Crab_850.fits", head_ja
mwrfits,stdmap, "ale2h_reso_xpol_Crab_850.fits", head_ja


end
