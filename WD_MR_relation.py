import glob as glob
import os
import sys
import numpy as np
import pandas as pd
from astropy import constants as const
import astropy.units as u
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import math

dirn=os.path.dirname(__file__)
sys.path.append(dirn)

#load Althaus 2013 phot. tables
MR_tab_a=pd.read_csv(dirn+"/MR_grids/"+"all_wds_logg_z02_althaus13.csv",sep=",")
MR_tab_ase=pd.read_csv(dirn+"/MR_grids/"+"all_wds_mass_z02_althaus13.csv",sep=",")

#load Bedard 2020 phot.constant mass tables   
MR_tab_m=pd.read_csv(dirn+"/MR_grids/"+"Bedard20_Table_mass_phot.csv",sep=",")
MR_tab_m["Lum"]=10**((MR_tab_m["Mbol"]-4.8)/(-2.5))
MR_tab_m["R"]=np.sqrt(MR_tab_m["Lum"])/(np.divide(MR_tab_m["Teff"],5770))**2

#load Bedard 2020 phot.constant logg tables   
MR_tab_s=pd.read_csv(dirn+"/MR_grids/"+"Bedard20_Table_spec_phot.csv",sep=",")
MR_tab_s["Lum"]=10**((MR_tab_s["Mbol"]-4.8)/(-2.5))
MR_tab_s["R"]=np.sqrt(MR_tab_s["Lum"])/(np.divide(MR_tab_s["Teff"],5770))**2

MR_tab_evol=pd.read_csv(dirn+"/MR_grids/"+"Bedard20_evol_seq.csv",sep=",")


def Radius_from_teff_logg(Temp,logg, Model, spec_type,thickness):
    """ 
    Input parameters 
        Teff(in K) and logg
        Model inputs: Althaus13 or Bedard20
        thickness: thick or thin 
        spec_type: DA or DB
    returns 
        radius and mass
    """
    Temp=np.array(Temp)
    logg=np.array(logg)
    Temp=np.log10(Temp)
    
    MR_tab_aseq=MR_tab_ase.loc[(MR_tab_ase["spec_type"]==spec_type)]
    MR_tab_an=MR_tab_a.loc[(MR_tab_a["spec_type"]==spec_type)]
    
    MR_tab_ms=MR_tab_m.loc[(MR_tab_m["spec_type"]==spec_type)]
    MR_tab_sp=MR_tab_s.loc[(MR_tab_s["spec_type"]==spec_type)]
    MR_tab_ev=MR_tab_evol.loc[(MR_tab_evol["thickness"]==thickness)]

    if Model=="Althaus13": 
        R=griddata((MR_tab_aseq["log(teff)"],MR_tab_aseq["logg"]),np.log10(MR_tab_aseq["R(Rsun)"]),(Temp,logg))
     
        R=np.where(np.isnan(R),griddata((MR_tab_an["log(teff)"],MR_tab_an["logg"]),np.log10(MR_tab_an["R(Rsun)"]),(Temp,logg)),R)
        
        if np.any(np.isnan(R)):
            print('-------------------------------')
            print('--> input parameters out of all the model grids for nan values')

    
    if Model=="Bedard20":
        R=griddata((MR_tab_ev["logT"],MR_tab_ev["logg"]),MR_tab_ev["logR"],(Temp,logg))

        R=np.where(np.isnan(R),griddata((np.log10(MR_tab_ms["Teff"]),MR_tab_ms["logg"]),np.log10(MR_tab_ms["R"]), (Temp,logg)),R)
             
        R=np.where(np.isnan(R),griddata((np.log10(MR_tab_sp["Teff"]),MR_tab_sp["logg"]),np.log10(MR_tab_sp["R"]),(Temp,logg)),R)
                
        if np.any(np.isnan(R)):
            print('-------------------------------')
            print('--> input parameters out of all the model grids for nan values')
    return 10**R

def Mass_from_teff_logg(Temp,logg, Model, spec_type,thickness):
    """ 
    Input parameters 
        Teff(in K) and logg
        Model inputs: Althaus13 or Bedard20
        thickness: thick or thin 
        spec_type: DA or DB
    returns 
        radius and mass
    """
    Temp=np.array(Temp)
    logg=np.array(logg)
    Temp=np.log10(Temp)
    
    MR_tab_aseq=MR_tab_ase.loc[(MR_tab_ase["spec_type"]==spec_type)]
    MR_tab_an=MR_tab_a.loc[(MR_tab_a["spec_type"]==spec_type)]
    
    MR_tab_ms=MR_tab_m.loc[(MR_tab_m["spec_type"]==spec_type)]
    MR_tab_sp=MR_tab_s.loc[(MR_tab_s["spec_type"]==spec_type)]
    MR_tab_ev=MR_tab_evol.loc[(MR_tab_evol["thickness"]==thickness)]

    if Model=="Althaus13": 
        M=griddata((MR_tab_aseq["log(teff)"],MR_tab_aseq["logg"]),MR_tab_aseq["mWD"],(Temp,logg))
       
        M=np.where(np.isnan(M),griddata((MR_tab_an["log(teff)"],MR_tab_an["logg"]),MR_tab_an["mWD"],(Temp,logg)),M)
        
    if Model=="Bedard20":
        M=griddata((MR_tab_ev["logT"],MR_tab_ev["logg"]),MR_tab_ev["Mass"],(Temp,logg))
       
        M=np.where(np.isnan(M),griddata((np.log10(MR_tab_ms["Teff"]),MR_tab_ms["logg"]),MR_tab_ms["Mass"], (Temp,logg)),M)
         
        M=np.where(np.isnan(M),griddata((np.log10(MR_tab_sp["Teff"]),MR_tab_sp["logg"]),MR_tab_sp["Mass"],(Temp,logg)),M)
        
        if np.any(np.isnan(M)):
            print('-------------------------------')
            print('--> input parameters out of all the model grids for nan values')
    return M
        
#example
Tempi=np.linspace(8000,70000,10)
loggi=np.linspace(7,9.3,10)
#print('Temp:'+str(Tempi))
#print('logg:'+str(loggi))

mon_r= Radius_from_teff_logg(Tempi, loggi,"Bedard20",'DA','thick')
#print(mon_r)
mon_m= Mass_from_teff_logg(Tempi, loggi,"Bedard20",'DA','thick')


arg_r= Radius_from_teff_logg(Tempi, loggi,"Althaus13",'DA','thick')
#print(arg_r)
arg_m= Mass_from_teff_logg(Tempi, loggi,"Althaus13",'DA','thick')



#plot DA models
plt.clf()
argen=MR_tab_ase.loc[(MR_tab_ase["spec_type"]=='DA')]
argen_s=MR_tab_a.loc[(MR_tab_a["spec_type"]=='DA')]
mont_g=MR_tab_m.loc[(MR_tab_m["spec_type"]=='DA')]
mont=MR_tab_s.loc[(MR_tab_s["spec_type"]=='DA')]
mont_ev=MR_tab_evol.loc[(MR_tab_evol["thickness"]=='thick')]

ax=plt.subplot(111)


#ax.plot(argen_s["logg"],10**argen_s["log(teff)"],'k.',alpha=0.5,label="Althaus DA logg")
#ax.plot(mont["logg"],mont["Teff"],'b.',label="Bedard DA logg")
ax.plot(mont_ev["logg"],10**mont_ev["logT"],'c.', label="Bedard evol. seq.")
ax.plot(mont_g["logg"],mont_g["Teff"],'y+',label="Bedard DA mass")
ax.plot(argen["logg"],10**argen["log(teff)"],'r.',label="Althaus DA mass")
# grp=argen.groupby(["spec_type","mWD"])
# grp1=mont_g.groupby(["spec_type","Mass"])
# for i in range(len(grp.first())):
#     ax.text(grp["logg"].first()[i]-0.07,10**grp["log(teff)"].first()[i]-400,grp["mWD"].first()[i],color='r',size=10)
# for i in range(len(grp1.first())):
#     ax.text(grp1["logg"].first()[i]-0.07,grp1["Teff"].first()[i]-280,grp1["Mass"].first()[i],color='k',size=10)
    
#ax.plot(loggi, Tempi, 'kx')
plt.xlabel("log g",size=14)
plt.ylabel("T$_{eff}$",size=14)
plt.yscale('log')
plt.legend(loc='upper right',ncol=1)
plt.show()
