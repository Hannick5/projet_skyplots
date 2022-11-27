import fonctions
import numpy as np
import matplotlib.pyplot as plt
import geopandas
from pathlib import Path
import lecture_de_TLE


# Définition des paramètres du satellite iss BIIF-9

a_iss = fonctions.calcul_a(0.0011264398169714)

i_iss = np.radians(51.6456)

e_iss = 0.0003349

w0 = np.radians(lecture_de_TLE.elements_keplerien("données/iss.txt")['Argument du périgée'][0])

W0 = np.radians(lecture_de_TLE.elements_keplerien("données/iss.txt")['Ascension'][0])

n0 = 0.0011264398169714

J2 = 1.08262652304819194e-3

M_iss = lecture_de_TLE.elements_keplerien("données/iss.txt")['Anomalie moyenne'][0]


vitesse_angulaire = -7.2921151467e-5

heure_tle = fonctions.tle_to_utc(float(lecture_de_TLE.elements_keplerien("données/iss.txt")['Heure'][0]))

date_tle = fonctions.t(fonctions.dateJulienne(2021,11,7,heure_tle))

DATE = np.linspace(fonctions.t(2459526.5),fonctions.t(2459527.5),num=10000)


X = []

Y = []


semaine_bis = ['08/11/2021','09/11/2021','10/11/2021','11/11/2021','12/11/2021','13/11/2021','14/11/2021']


semaine = [2459526.5,2459527.5,2459528.5,2459529.5,2459530.5,2459531.5,2459533.49931]

for k in range(1,len(semaine)):
    DATE = np.linspace(fonctions.t(semaine[k-1]),fonctions.t(semaine[k]),num=10000)
    for date in DATE:
        w_iss = fonctions.var_omega(w0,n0,J2,a_iss,e_iss,i_iss,date,date_tle)

        W_iss = fonctions.var_OMEGA(W0,n0,J2,a_iss,e_iss,i_iss,date,date_tle)

        M_corr_iss = fonctions.var_anomalie_moyenne(M_iss,n0,J2,a_iss,e_iss,i_iss,date,date_tle)

        E_iss = fonctions.kepler_resolve(M_corr_iss,e_iss,1e-6)

        coord_orbital_iss = fonctions.coord_orbital(a_iss,e_iss,E_iss)

        coord_celeste_iss = fonctions.coord_celeste(i_iss,w_iss,W_iss,e_iss,E_iss,a_iss)

        coord_cart_iss = fonctions.coord_cart(i_iss,w_iss,W_iss,e_iss,E_iss,a_iss,date,vitesse_angulaire)

        coord_geodesique_iss = fonctions.coord_geodesique(i_iss,w_iss,W_iss,E_iss,date,vitesse_angulaire)

        X.append(np.degrees(coord_geodesique_iss[0]))
        Y.append(np.degrees(coord_geodesique_iss[1]))



    world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))

    fig, ax = plt.subplots()
    world.plot(ax=ax)
    plt.scatter(X,Y,s=10,c='red')
    plt.title(f"Trace de l'ISS de {semaine_bis[0]} à {semaine_bis[k]}")
    plt.xlabel('Longitude(deg)')
    plt.ylabel('Latitude(deg)')
    plt.show(block=False)
    plt.pause(3)
    plt.close()

