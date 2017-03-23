import matplotlib.pyplot as plt
import numpy as np
# Results

# 

qsolist = []
qsolist.append({'QSO':'3C286','EXP':'XPOL','FREQ':'86 GHz',  'p':13.5,'ep':0.3,'ang':37.3,'eang':0.8,'color':'hotpink'})
qsolist.append({'QSO':'3C286','EXP':'NIKA','FREQ':'150 GHz', 'p':13.6,'ep':0.8,'ang':28.0,'eang':0.9+1.8,'color':'hotpink'})
qsolist.append({'QSO':'3C286','EXP':'CARMA','FREQ':'225 GHz','p':0.0,'ep':0.0,'ang':39.1,'eang':1.0,'color':'hotpink'})
qsolist.append({'QSO':'3C286','EXP':'XPOL','FREQ':'229 GHz', 'p':14.4,'ep':1.8,'ang':33.1,'eang':5.7,'color':'hotpink'})
qsolist.append({'QSO':'3C286','EXP':'SMA 2006','FREQ':'230 GHz','p':11.5,'ep':2.5,'ang':35.6,'eang':5.9,'color':'hotpink'})
qsolist.append({'QSO':'3C286','EXP':'ALMA','FREQ':'230 GHz', 'p':16.7,'ep':0.2,'ang':39,'eang':0.4,'color':'hotpink'})
qsolist.append({'QSO':'3C286','EXP':'NIKA','FREQ':'260 GHz', 'p':14.3,'ep':1.7,'ang':30.0,'eang':2.5+1.8,'color':'hotpink'})
qsolist.append({'QSO':'3C286','EXP':'SMA 2016','FREQ':'340 GHz','p':15.7,'ep':0.8,'ang':37.4,'eang':1.5,'color':'hotpink'})
qsolist.append({'QSO':'3C273','EXP':'XPOL','FREQ':'86 GHz', 'p':1.1,'ep':0.1,'ang':-37.8,'eang':0.9,'color':'deepskyblue'})
qsolist.append({'QSO':'3C273','EXP':'NIKA','FREQ':'150 GHz','p':2.0,'ep':0.1,'ang':-74.05,'eang':1.1+1.8,'color':'deepskyblue'})
qsolist.append({'QSO':'3C273','EXP':'XPOL','FREQ':'229 GHz','p':3.6,'ep':0.2,'ang':-76.8,'eang':1.6,'color':'deepskyblue'})
qsolist.append({'QSO':'3C273','EXP':'NIKA','FREQ':'260 GHz','p':3.41,'ep':0.24,'ang':-88.7,'eang':0.9+1.8,'color':'deepskyblue'})
qsolist.append({'QSO':'3C273','EXP':'XPOL','FREQ':'260 GHz','p':1.59,'ep':0.18,'ang':-71.8,'eang':3.1,'color':'deepskyblue'})
qsolist.append({'QSO':'3C279','EXP':'NIKA','FREQ':'150 GHz','p':9.5,'ep':0.6,'ang':31.9,'eang':0.7+1.8,'color':'limegreen'})
qsolist.append({'QSO':'3C279','EXP':'NIKA','FREQ':'260 GHz','p':9.8,'ep':0.4,'ang':35.88,'eang':0.5+1.8,'color':'limegreen'})
qsolist.append({'QSO':'3C279','EXP':'XPOL','FREQ':'260 GHz','p':11.79,'ep':0.29,'ang':45.6,'eang':0.7,'color':'limegreen'})
qsolist.append({'QSO':'0923+392','EXP':'NIKA','FREQ':'150 GHz','p':2.72,'ep':0.24,'ang':-50.36,'eang':1.74+1.8,'color':'violet'})
qsolist.append({'QSO':'0923+392','EXP':'NIKA','FREQ':'260 GHz','p':3.25,'ep':0.28,'ang':-46.1,'eang':2.4+1.8,'color':'violet'})
qsolist.append({'QSO':'0923+392','EXP':'XPOL','FREQ':'260 GHz','p':6.1,'ep':2.3,'ang':-52.6,'eang':10.0,'color':'violet'})

nlist = len(qsolist)
dp =   []
dang = []
dnp = []
fr  = []
dnang = []
dcp=[]
dcang=[]

for qso in qsolist:
    name = qso['EXP']
    freq = qso['FREQ']
    me = qso['p']
    ms = qso['ep']
    color = qso['color']
    if ms > 0.0:
        dp.append(np.random.normal(me,ms,100))
        dnp.append(name)
        fr.append(freq)
        dcp.append(color)
        me = qso['ang']
        ms = qso['eang']
    if ms > 0.0:
        dang.append(np.random.normal(me,ms,1000))
        dnang.append(name)
        dcang.append(color)
        
qsoname=['3C286','3C273','3C279','0923+392']

fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(16, 12))
bplot=axes.boxplot(dp,   showfliers=False,
                         vert=True,   # vertical box aligmnent
                         patch_artist=True)
for patch, color in zip(bplot['boxes'], dcp):
    patch.set_facecolor(color)

plt.figtext(0.2,0.8,qsoname[0], color='hotpink',fontsize=20,fontweight='bold')
plt.figtext(0.45,0.3,qsoname[1],color='deepskyblue',fontsize=20,fontweight='bold')
plt.figtext(0.65,0.65,qsoname[2],color='limegreen',fontsize=20,fontweight='bold')
plt.figtext(0.73,0.4,qsoname[3],color='violet',fontsize=20,fontweight='bold')

plt.grid()
plt.ylabel('Polarization degree [%]',fontsize=24)
axes.set_ylabel('Polarization degree [%]',fontsize=24)
axes.set_xlim(0,19)
axes.tick_params(labelsize=18)
plt.setp(axes, xticks=[y+1 for y in range(len(fr))],
         xticklabels=fr)
axes.xaxis.set_ticks_position('bottom')
plt.xticks(rotation='vertical')

axes2 = axes.twiny()
axes2.set_xlim(0,19)
plt.setp(axes2,xticks=[y+1 for y in range(len(fr))],xticklabels=dnp)
axes2.xaxis.set_ticks_position('top')
axes2.tick_params(labelsize=16.5)
plt.xticks(rotation='vertical')

fig2, axes3 = plt.subplots(nrows=1, ncols=1, figsize=(16, 12))
bplot2=axes3.boxplot(dang,   showfliers=False,# notch shape
                         vert=True,   # vertical box aligmnent
                         patch_artist=True)

for patch, color in zip(bplot2['boxes'], dcang):
    print color
    patch.set_facecolor(color)

plt.setp(axes3, xticks=[y+1 for y in range(len(fr))],
         xticklabels=fr)
axes3.set_xlim(0,19)
axes3.xaxis.set_ticks_position('bottom')
axes3.tick_params(labelsize=18)
plt.xticks(rotation='vertical')

plt.figtext(0.2,0.65,qsoname[0], color='hotpink',fontsize=20,fontweight='bold')
plt.figtext(0.45,0.35,qsoname[1],color='deepskyblue',fontsize=20, fontweight='bold')
plt.figtext(0.65,0.65,qsoname[2],color='limegreen',fontsize=20,fontweight='bold')
plt.figtext(0.73,0.50,qsoname[3],color='violet',fontsize=20,fontweight='bold')

plt.grid()
plt.ylabel('Polarization angle [degrees]',fontsize=24)
axes3.set_ylabel('Polarization angle [degrees]',fontsize=24)
axes4 = axes3.twiny()
axes4.set_xlim(0,19)
plt.setp(axes4,xticks=[y+1 for y in range(len(fr))],xticklabels=dnang)
#plt.setp(xticknames, fontsize=20)
axes4.xaxis.set_ticks_position('top')
axes4.tick_params(labelsize=16.5)
plt.xticks(rotation='vertical')

plt.show()







