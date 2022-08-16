# Bibliotheken Import
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import pandas as pd

# Funktionsdefinition
cos  = lambda arg : np.cos(np.deg2rad(arg))
sin  = lambda arg : np.sin(np.deg2rad(arg))
acos = lambda arg : np.rad2deg(np.arccos(arg))
asin = lambda arg : np.rad2deg(np.arcsin(arg))

# Parameterdefinition
bg = 47.1733                        # Breitengrad
lg = 9.4727                         # Längengrad
al = 0.3                            # Albedo-Faktor
pr = 0.8                            # Performance Ration 
bf = 0.9                            # Bifacial-Faktor
aziFlArray = [0, 90, 180, -90]      # [grad] Süden, Westen, Norden, Osten
neigArray = np.linspace(0, 90, 90)  # [grad] Austellungswinkel

#df = pd.read_csv("data.csv", delimiter=';', header = 0)
df = pd.read_csv("data.csv", delimiter=';', header = 0)
df = df.set_index(pd.DatetimeIndex(df['Time'])) 

hGlo  = df['hGlo'].values   # [W/m2]
hDif  = df['hDif'].values   # [W/m2]
hDNI = hGlo - hDif          # [W/m2]
tAmb  = df['Tamb'].values   # [°C]
tutc = df.index             # [W/m2]

# Delta Zeit berechnen
lfStd = np.zeros(tutc.size)
for t in range(tutc.size):
    # berechnet laufender Tag im Jahr
    noDay = (tutc[t] - dt.datetime(tutc[0].year, 1, 1, 0)).days 
    # [h] berechnet laufende Stunde im Tag
    noHou = tutc[t].hour + (tutc[t].minute)/60 + (tutc[t].second)/3600 
    lfStd[t] = noDay*24 + noHou
deltaT = lfStd[1] - lfStd[0] # [h]

# Sonnenstand
omega = 15 * np.mod(lfStd, 24) + lg - 180 # [grad] Stundenwinkel
dekl = 23.45 * cos(360 / 8760 * (lfStd - 173 * 24)) #[grd] Dekliation
h = asin(sin(dekl) * sin(bg) + cos(dekl) * cos(bg) * cos(omega))
azi = acos((sin(h) * sin(bg) - sin(dekl)) / (cos(h) * cos(bg))) * np.sign(omega)

################# Berechnung Verschiedene Ausrichtungen #############################

for aziFl in aziFlArray:

    ################# Normales Photovoltaikmodule #############################

    # JahresStrahlungsenergie
    hw = np.zeros(neigArray.size)
    for n in range(neigArray.size):
        neig = neigArray[n]

        # 3_k-Model
        cosTheta = cos(neig) * sin(h) + sin(neig) * cos(h) * cos(azi - aziFl)
        cosTheta[cosTheta < 0] = 0 # alle winkel über neunzig Grad werden null gesetzt
        hDirFl = hDNI * cosTheta
        hDifFl = hDif * 0.5 * (1 + cos(neig))
        hAlbFl = hGlo * al * 0.5 * (1 - cos(neig))
        hFl = hDirFl + hDifFl + hAlbFl # [W/m^2]
        hw[n] = pr * np.sum(hFl)*deltaT/1000 # [kWh/kW] pro Killowat peal elektrische Energie

    # Darstellung
    plt.plot(neigArray, hw , 'b-', label ='PV')

    ################# Bifaciale Photovoltaikmodule #############################

    # JahresStrahlungsenergie berechnen
    hwb = np.zeros(neigArray.size)
    for n in range(neigArray.size):
        neig = neigArray[n]

        # 6_k-Model
        #oben
        cosThetab = cos(neig) * sin(h) + sin(neig) * cos(h) * cos(azi - aziFl)
        cosThetab[cosThetab < 0] = 0 #alle winkel über neunzig Grad werden null gesetzt
        hDirFlbo = hDNI * cosThetab
        hDifFlbo = hDif * 0.5 * (1 + cos(neig))
        hAlbFlbo = hGlo * al * 0.5 * (1 - cos(neig))
        #unten
        cosThetabu = - (cos(neig) * sin(h) + sin(neig) * cos(h) * cos(azi - aziFl)) # entweder neigung um 90* drehen oder forzeichen Tauschen wegen dem aufstrahlwingel
        cosThetabu[cosThetabu < 0] = 0 #alle winkel über neunzig Grad werden null gesetzt
        hDirFlbu = hDNI * cosThetabu * bf
        hDifFlbu = hDif * 0.5 * (1 + cos(neig)) * bf
        hAlbFlbu = hGlo * al * 0.5 * (1 - cos(neig)) * bf

        hFlb = hDirFlbo + hDifFlbo + hAlbFlbo + hDirFlbu + hDifFlbu + hAlbFlbu # [W/m^2]
        hwb[n] = pr * np.sum(hFlb)*deltaT/1000 # [kWh/kW] pro Killowat peal elektrische Energie

    # Darstellung
    plt.plot(neigArray, hwb , 'g-', label ='Bifaciale PV')
    plt.xlabel('Neigung [grad]')
    plt.ylabel('Jahresstrahlungsenergie [kWh/kW]')
    plt.legend(loc='best')
    plt.title(f'Ausrichtung {aziFl}[Grad]')
    plt.grid()
    plt.show()
