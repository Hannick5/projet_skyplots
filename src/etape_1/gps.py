import fonctions
import numpy as np
import matplotlib.pyplot as plt
import geopandas
from pathlib import Path
import lecture_de_TLE

# Définition des paramètres du satellite GPS BIIF-9

a_gps = fonctions.calcul_a(1.45*(10**(-4)))

i_gps = np.radians(53.8532)

e_gps = 0.0063356

n0 = 1.45*(10**(-4))

J2 = 1.08262652304819194e-3

w0 = np.radians(lecture_de_TLE.elements_keplerien("données/gps-ops.txt")['Argument du périgée'][22])

W0 = np.radians(lecture_de_TLE.elements_keplerien("données/gps-ops.txt")['Ascension'][22])

M_gps = lecture_de_TLE.elements_keplerien("données/gps-ops.txt")['Anomalie moyenne'][22]

vitesse_angulaire = -7.2921151467e-5

heure_tle = fonctions.tle_to_utc(float(lecture_de_TLE.elements_keplerien("données/gps-ops.txt")['Heure'][22]))

date_tle = fonctions.t(fonctions.dateJulienne(2021,11,8,heure_tle))

DATE = np.linspace(fonctions.t(2459526.5),fonctions.t(2459527.5),num=10000)

X = []

Y = []

semaine_bis = ['08/11/2021','09/11/2021','10/11/2021','11/11/2021','12/11/2021','13/11/2021','14/11/2021']


semaine = [2459526.5,2459527.5,2459528.5,2459529.5,2459530.5,2459531.5,2459533.49931]

for k in range(1,len(semaine)):
    DATE = np.linspace(fonctions.t(semaine[k-1]),fonctions.t(semaine[k]),num=10000)
    for date in DATE:
        w_gps = fonctions.var_omega(w0,n0,J2,a_gps,e_gps,i_gps,date,date_tle)

        W_gps = fonctions.var_OMEGA(W0,n0,J2,a_gps,e_gps,i_gps,date,date_tle)

        M_corr_gps = fonctions.var_anomalie_moyenne(M_gps,n0,J2,a_gps,e_gps,i_gps,date,date_tle)

        E_gps = fonctions.kepler_resolve(M_corr_gps,e_gps,1e-6)

        coord_orbital_gps = fonctions.coord_orbital(a_gps,e_gps,E_gps)

        coord_celeste_gps = fonctions.coord_celeste(i_gps,w_gps,W_gps,e_gps,E_gps,a_gps)

        coord_cart_gps = fonctions.coord_cart(i_gps,w_gps,W_gps,e_gps,E_gps,a_gps,date,vitesse_angulaire)

        coord_geodesique_gps = fonctions.coord_geodesique(i_gps,w_gps,W_gps,E_gps,date,vitesse_angulaire)

        X.append(np.degrees(coord_geodesique_gps[0]))
        Y.append(np.degrees(coord_geodesique_gps[1]))


    world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))

    fig, ax = plt.subplots()
    world.plot(ax=ax)
    plt.scatter(X,Y,s=10,c='red')
    plt.title(f'Trace du satellite GPS BIIF-9 de {semaine_bis[0]} à {semaine_bis[k]}')
    plt.xlabel('Longitude(deg)')
    plt.ylabel('Latitude(deg)')
    plt.show(block=False)
    plt.pause(3)
    plt.close()

