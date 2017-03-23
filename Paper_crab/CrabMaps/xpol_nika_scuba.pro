;; XPOL

dir_xpol = '/Users/ritacco/Desktop/maps/'

map_ja_i = mrdfits(dir_xpol+'crab_I_coadd28_reproj_nbl_ccentonly.fits',0,head_ja)
map_ja_q = mrdfits(dir_xpol+'crab_Q_coadd28_reproj_nbl_ccentonly.fits',0,head_ja)
map_ja_u = mrdfits(dir_xpol+'crab_U_coadd28_reproj_nbl_ccentonly.fits',0,head_ja)
extast,head_ja,astr_ja
astr_ja.ctype[0] = "RA---TAN"
astr_ja.ctype[1] = "DEC--TAN"

reso = astr_ja.cdelt[1]*3600.0

longobj = sxpar(head_ja,'crval1') 
latobj = sxpar(head_ja,'crval2') 
ra_xpol = SIXTY(longobj/15.0)
dec_xpol = SIXTY(latobj)
dec_xpol[2] = Float(Round(dec_xpol[2]*1000)/1000.)
message, /info, 'IMBFITS pointing is :'
message, /info, 'R.A.: '+strtrim(ra_xpol[0],2)+' h '+strtrim(ra_xpol[1],2)+' min '+strtrim(ra_xpol[2],2)+' s'
message, /info, 'Dec.: '+strtrim(dec_xpol[0],2)+' deg '+strtrim(dec_xpol[1],2)+" arcmin "+strtrim(dec_xpol[2],2)+' arcsec'

coord_plot_xpol = [ten(ra_xpol[0],ra_xpol[1],ra_xpol[2])*15.0,$
                    ten(dec_xpol[0],dec_xpol[1],dec_xpol[2])]    

ra_xpol  = sxpar(head_ja, 'crval1')
dec_xpol = sxpar(head_ja, 'crval2')
nx = sxpar(head_ja, 'naxis1')
ny = sxpar(head_ja, 'naxis2')
reso = sxpar(head_ja, "CDELT2", /silent) * 3600.d0
;; reso = abs(sxpar(head_ja, "CD1_1", /silent)) * 3600.d0
naxis = [nx, ny]                                       ;Nb pixel along x and y
cd = [[1.0,-0.0],[0.0,1.0]]                            ;Rotation matrix but no rotation here
cdelt = [-1.0, 1.0] * reso/3600.0                      ;Pixel size (ra along -1)
crpix = [37.0,20.0]                                    ;Ref pixel (central pixel (always odd nb))
ra = ra_xpol                                          ;RA in degrees
dec = dec_xpol                                        ;DEC in degrees
crval = [ra_xpol, dec_xpol]                                      ;ra dec of the ref pix
ctype = ['RA---TAN','DEC--TAN']                        ;Projection type
     
astrometry = {naxis:naxis, $
                   cd:cd, $
                   cdelt:cdelt, $
                   crpix:crpix, $
                   crval:crval, $
                   ctype:ctype,$
                   longpole:180.0, $
                   latpole:90.0, $
                   pv2:[0.0,0.0]}


mkhdr,  head_xpol, map_ja_i ;get header typique
putast, head_xpol, astrometry, equinox=2000, cd_type=0 ;astrometry in header

map_ja_i = map_ja_i*6.; Jy factor conversion Aumont et al.
map_ja_q = map_ja_q*6.; Jy factor conversion Aumont et al.
map_ja_u = map_ja_u*6.; Jy factor conversion Aumont et al.

dir = '/Users/ritacco/Nika/software/Processing/Papers_Proposal_Notes/Paper_crab/CrabMaps/'

nika_ar        = mrdfits(dir+'ale2h_reso_xpol_Crab_2mm.fits', 1, header)
nika_ar_stddev = mrdfits(dir+'ale2h_reso_xpol_Crab_2mm.fits', 2, head)
map_i_2mm = nika_ar[*,*,0]
map_q_2mm = nika_ar[*,*,1]
map_u_2mm = nika_ar[*,*,2]

map_i_2mm = nika_ar[*,*,0]
map_q_2mm = nika_ar[*,*,1]
map_u_2mm = nika_ar[*,*,2]

map_i_2mm_stddev = nika_ar_stddev[*,*,0]
map_q_2mm_stddev = nika_ar_stddev[*,*,1]
map_u_2mm_stddev = nika_ar_stddev[*,*,2]

extast,header,astr
omega_b_2mm = 2*!pi*(18.2*!fwhm2sigma)^2
ipol_2mm  = sqrt(map_q_2mm^2+map_u_2mm^2)
ipol_xpol = sqrt(map_ja_q^2+map_ja_u^2)

i_2mm  = filter_image(map_i_2mm, fwhm=sqrt(25.2^2-18.2^2)/reso)*1.28*reso*reso/omega_b_2mm
ip_2mm = filter_image(ipol_2mm,  fwhm=sqrt(25.2^2-18.2^2)/reso)*1.28*reso*reso/omega_b_2mm

spec_i   = [alog(i_2mm/map_ja_i)/alog(150./90.)]
spec_ip  = [alog(ip_2mm/ipol_xpol)/alog(150./90.)]


head = head_ja
sxaddpar, head, "naxis", 2
sxdelpar, head, "naxis3"
sxaddpar, head, "ctype1", "RA---TAN"
sxaddpar, head, "ctype2", "DEC--TAN"

levels=[0.2, 0.6, 1, 1.4, 1.8, 2.2]

ar_my_imview, map_i_2mm, header=head_xpol, contour=filter_image(map_ja_i,fwhm=3), c_header=head_xpol, $
              levels=levels






end
