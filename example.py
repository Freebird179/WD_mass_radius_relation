import WD_MR_relation as MR
import matplotlib.pyplot as plt

#example
Tempi=[8560,12500,28300]
loggi=[7.35,9.21,8.56]
print('Temp:'+str(Tempi))
print('logg:'+str(loggi))

mon_r= MR.Radius_from_teff_logg(Tempi, loggi,"Bedard20",'DA','thick')
print('Radius_mont :'+str(mon_r))
mon_m= MR.Mass_from_teff_logg(Tempi, loggi,"Bedard20",'DA','thick')
print('Mass_mont :'+str(mon_m))

arg_r= MR.Radius_from_teff_logg(Tempi, loggi,"Althaus13",'DA','thick')
print('Radius_arg:'+str(arg_r))
arg_m= MR.Mass_from_teff_logg(Tempi, loggi,"Althaus13",'DA','thick')
print('Mass_arg:'+str(arg_m))
