import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
import warnings
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import math as m
from astropy import units as u 
from specutils.fitting import find_lines_derivative
from specutils import Spectrum1D
from specutils import SpectralRegion
from specutils.fitting import fit_generic_continuum
from specutils.manipulation import extract_region
from pandas import read_csv
from plotly.subplots import make_subplots
from astropy.coordinates import get_constellation, SkyCoord

#On posa Espectres cal que es posi la direcció on es troba l'arxiu mastar-goodspec-v3_1_1-v1_7_7.fits entre les parèntesis
#On posa Gaia_Data cal que es posi la direcció on es troba l'arxiu mastarall-gaiadr2-extcorr-simbad-ps1-v3_1_1-v1_7_7-v1.fits entre les parèntesis
Espectres = fits.open('C:\\Users\\Usuari\\Desktop\\Tdr\\mastar-goodspec-v3_1_1-v1_7_7.fits')[1].data
Gaia_Data = fits.open('C:\\Users\\Usuari\\Downloads\\mastarall-gaiadr2-extcorr-simbad-ps1-v3_1_1-v1_7_7-v1.fits')[1].data
id = str(input("Introdueix MANGAID: "))
star_id = (np.where(Espectres["MANGAID"] == id)[0])[0]
(np.where(Espectres["MANGAID"] == id)[0])[0]
Longitud_Ona = Espectres["WAVE"][star_id] *u.AA
Flux_Estrella = Espectres["Flux"][star_id]
Gaia_index = int((np.where(Gaia_Data["PSFMAG"] == Espectres["PSFMAG"][star_id])[0])[0])
Gaia_Data_Py = pd.DataFrame(Gaia_Data['R_EST'])
Distancia = int(Gaia_Data_Py.iloc[Gaia_index])

if Espectres["PHOTOCAT"][star_id] == "sdss_dr8":
  EspectresUGRIZ =  [float(i) for i in Espectres["PSFMAG"][star_id]]
  g = EspectresUGRIZ[1]
  r = EspectresUGRIZ[2]
  B = float(g + 0.3130*(g - r) + 0.2271)
  V = float(g - 0.5784*(g - r) - 0.0038)
  Color_Index = B - V
elif Espectres["PHOTOCAT"][star_id] == "ps1":
    EspectresGRIZ = [float(i) for i in Espectres["PSFMAG"][star_id]]
    g = EspectresGRIZ[1]
    r = EspectresGRIZ[2]
    B = float(g + 0.3130*(g - r) + 0.2271)
    V = float( g - 0.5784*(g - r) - 0.0038)
    Color_Index = B - V
elif Espectres["PHOTOCAT"][star_id] == "apass_dr8":
    EspectresGRI = [float(i) for i in Espectres["PSFMAG"][star_id]]
    g = EspectresGRI[1]
    r = EspectresGRI[2]
    B = float(g + 0.3130*(g - r) + 0.2271)
    V = float(g - 0.5784*(g - r) - 0.0038)
    Color_Index = B - V
elif Espectres["PHOTOCAT"][star_id] == "gaia_dr2":
    EspectresGAIA = [float(i) for i in Espectres["PSFMAG"][star_id]]
    G = float(EspectresGAIA[1])
    G_BP = float(EspectresGAIA[2])
    G_RP = float(EspectresGAIA[3])
    V = 0.02704 - 0.01424*(G_BP - G_RP) + 0.2156*((G_BP - G_RP)**2.0) - 0.01426*((G_BP - G_RP)**3.0) + G
    g = - 0.2199 + 0.6365*(G_BP - G_RP) + 0.1548*((G_BP - G_RP)**2.0) - 0.0064*((G_BP - G_RP)**3.0) + G
    r = +0.09837 - 0.08592*(G_BP - G_RP) - 0.1907*((G_BP - G_RP)**2.0) + 0.1701*((G_BP - G_RP)**3.0) - 0.02263*((G_BP - G_RP)**4.0) + G
    B = g + 0.3130*(g - r) + 0.2271
    Color_Index = B - V

Teff = int(((4600*((1/(0.92 * Color_Index + 1.7)) + (1/(0.92 * Color_Index + 0.62))))))
Log_teff = m.log10(Teff)

BC = float(-8.499*((Log_teff-4)**4) + 13.421*((Log_teff-4)**3) - 8.131*((Log_teff-4)**2) - 3.901*(Log_teff-4) - 0.438)
MV = V - 2.5*((m.log10(Distancia/10))**2)
MBol = MV + BC

Gaia_Data_Pd = pd.DataFrame(Gaia_Data['R_EST'])
Distancia = int(Gaia_Data_Py.iloc[Gaia_index])

L = ((10**(-0.4*(MBol)))*((3.0128*(10**28))))
L_Sol = L/(3.827*(10**26))

Radi =  m.sqrt((L/(4*m.pi*(5.6703*(10**-8))*(Teff**4))))/(10**3)

Massa = L_Sol**(1/3.5)

Ra = Espectres["RA"][star_id]
Dec = Espectres["DEC"][star_id]
Coordinades = SkyCoord(Ra, Dec, unit = u.deg)
Constelació = get_constellation(Coordinades, short_name=True, constellation_list='iau')

Velocitat_Heliocentrica = Espectres["HELIOV"][star_id]

Classes_estel·lars = pd.read_csv('C:\\Users\\Usuari\\Downloads\\Stellar_classes.csv') #Entre les parèntesis cal que es posi la direcció de la taula Stellar_classes
Classes_estel·lar_1 = Classes_estel·lars.iloc[(Classes_estel·lars['Temp']-Teff).abs().argsort()[:3]]
Classe_Estrella = Classes_estel·lar_1.iloc[(Classes_estel·lar_1['Abs_Mag']-MV).abs().argsort()[:1]]['Stellar_Type']

EspectreD1 = Spectrum1D(Espectres["FLUX"][star_id]*(u.Jy), Espectres["WAVE"][star_id]*u.AA)
EspectreD1_net = extract_region(EspectreD1, SpectralRegion(4000* u.AA, 7000 * u.AA))
EspectreD1_X = EspectreD1_net.spectral_axis
EspectreD1_Y = EspectreD1_net.flux
Espectre_Ajustat = fit_generic_continuum(EspectreD1_net)
Flux_Ajustat = Espectre_Ajustat(EspectreD1_X)
Espectre_normalitzat = EspectreD1_net / Flux_Ajustat

n = -15
with warnings.catch_warnings():  # Ignora avisos
    warnings.simplefilter('ignore')
    linies = find_lines_derivative(Espectre_normalitzat, flux_threshold= n)
while len(linies[linies['line_type'] == 'absorption']) >= 30:
    n += 0.1
    with warnings.catch_warnings():
         warnings.simplefilter('ignore')
         linies = find_lines_derivative(Espectre_normalitzat, flux_threshold= n)
LiniesÀtoms = read_csv('C:\\Users\\Usuari\\Desktop\\Tdr\\ElementsTDR.csv') # Entre les parèntesis cal que es posi la direcció de la taula ElementsTDR
LiniesEspectre = linies[linies['line_type'] == 'absorption']['line_center'] / u.AA

LiniesNom = []
LiniesTrobades = []
Id_estrlla = str(Espectres["PLATE"][star_id])
for linies_abs in LiniesEspectre:
    for Linia_atom in LiniesÀtoms['Rest λ (Å)']:
            x = float(Linia_atom) - float(linies_abs)
            if  x <= 1.5 and x >= -1.5:
                n = LiniesÀtoms[LiniesÀtoms['Rest λ (Å)'] == Linia_atom].index[0]
                LiniesNom.append(LiniesÀtoms['Name'][n])
                LiniesTrobades.append(LiniesÀtoms['Rest λ (Å)'][n])

LiniesTrobades_Python = np.array(LiniesTrobades)


fig = px.line(x=EspectreD1_net.spectral_axis, y=EspectreD1_net.flux, title="Espectre estrella")
fig.update_yaxes(title_text='Fluxe (1e-17 erg/s/cm^2/Angstrom)')
fig.update_xaxes(title_text='Longitud Ona (Å)')

for LiniaTrobada in LiniesTrobades_Python:
        Valor_index = LiniesTrobades_Python.tolist().index(LiniaTrobada)
        LiniaNom = LiniesNom[Valor_index]
        Posició_Nom = int(LiniaTrobada-2)
        fig.add_vline(x=LiniaTrobada, line_width=1, line_dash="dash", line_color="orange", annotation_text=LiniaNom, annotation_position="bottom right")


fig1 = go.Figure(data=[go.Table(header=dict(values=["Classe estel·lar","Magnitud aparent", "Magnitud Absoluta", "Magnitud Bolometrica","Temperatura Efectiva", "Lluminositat(☉)", "Distancia(parsec)", "Radi(Km)","Massa(☉)", "Ra","Dec","Constel·lació","Velocitat Heliocèntrica(Km/s)"],
font=dict(size=10), align="left"), cells=dict(values=[Classe_Estrella,round(V,2), round(MV,2), round(MBol,2),round(Teff, 2), round(L_Sol,2), round(Distancia,2), round(Radi,2), round(Massa,2), round(Ra,2), Dec, Constelació, round(Velocitat_Heliocentrica, 2)], align = "left"))])


fig.update_layout(
    height=500,
    showlegend=False,
    title_text=("Espectre estrella: " + Espectres['MANGAID'][star_id]),
)
fig1.update_layout(
    height=500,
    showlegend=False,
    title_text=("Característiques estrella: " + Espectres['MANGAID'][star_id]),
)
fig.show()
fig1.show()
