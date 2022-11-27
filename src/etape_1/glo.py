import fonctions
import numpy as np
import matplotlib.pyplot as plt
import geopandas
from pathlib import Path
import lecture_de_TLE


# Définition des paramètres du satellite glo BIIF-9

a_glo = fonctions.calcul_a(1.45*(10**(-4)))

i_glo = np.radians(65.9702)

e_glo = 0.0025200

J2 = 1.08262652304819194e-3

w0 = np.radians(lecture_de_TLE.elements_keplerien("données/glo-ops.txt")['Argument du périgée'][8])

W0 = np.radians(lecture_de_TLE.elements_keplerien("données/glo-ops.txt")['Ascension'][8])

M_glo = lecture_de_TLE.elements_keplerien("données/glo-ops.txt")['Anomalie moyenne'][8]

n0 = 1.45*(10**(-4))

vitesse_angulaire = -7.2921151467e-5

heure_tle = fonctions.tle_to_utc(float(lecture_de_TLE.elements_keplerien("données/glo-ops.txt")['Heure'][8]))

date_tle = fonctions.t(fonctions.dateJulienne(2021,11,7,heure_tle))


DATE = np.linspace(fonctions.t(2459526.5),fonctions.t(2459527.5),num=10000)

X = []

Y = []



semaine_bis = ['08/11/2021','09/11/2021','10/11/2021','11/11/2021','12/11/2021','13/11/2021','14/11/2021']


semaine = [2459526.5,2459527.5,2459528.5,2459529.5,2459530.5,2459531.5,2459533.49931]

for k in range(1,len(semaine)):
    DATE = np.linspace(fonctions.t(semaine[k-1]),fonctions.t(semaine[k]),num=10000)
    for date in DATE:
        w_glo = fonctions.var_omega(w0,n0,J2,a_glo,e_glo,i_glo,date,date_tle)

        W_glo = fonctions.var_OMEGA(W0,n0,J2,a_glo,e_glo,i_glo,date,date_tle)

        M_corr_glo = fonctions.var_anomalie_moyenne(M_glo,n0,J2,a_glo,e_glo,i_glo,date,date_tle)

        E_glo = fonctions.kepler_resolve(M_corr_glo,e_glo,1e-6)

        coord_orbital_glo = fonctions.coord_orbital(a_glo,e_glo,E_glo)

        coord_celeste_glo = fonctions.coord_celeste(i_glo,w_glo,W_glo,e_glo,E_glo,a_glo)

        coord_cart_glo = fonctions.coord_cart(i_glo,w_glo,W_glo,e_glo,E_glo,a_glo,date,vitesse_angulaire)

        coord_geodesique_glo = fonctions.coord_geodesique(i_glo,w_glo,W_glo,E_glo,date,vitesse_angulaire)

        X.append(np.degrees(coord_geodesique_glo[0]))
        Y.append(np.degrees(coord_geodesique_glo[1]))

    world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))

    fig, ax = plt.subplots()
    world.plot(ax=ax)
    plt.scatter(X,Y,s=10,c='red')
    plt.title(f'Trace du satellite Glonass COSMOS 2459 du {semaine_bis[0]} à {semaine_bis[k]}')
    plt.xlabel('Longitude(deg)')
    plt.ylabel('Latitude(deg)')
    plt.show(block=False)
    plt.pause(3)
    plt.close()

