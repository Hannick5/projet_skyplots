import fonctions
import numpy as np
import matplotlib.pyplot as plt
import geopandas
import lecture_de_TLE

# Définition des paramètres du satellite galileo GSAT0209

a_galileo = fonctions.calcul_a(1.45*(10**(-4)))

i_galileo = np.radians(55.0919)

e_galileo = 0.0002373

n0 = 1.45*(10**(-4))

J2 = 1.08262652304819194e-3

w0 = np.radians(lecture_de_TLE.elements_keplerien("données/galileo.txt")['Argument du périgée'][11])

W0 = np.radians(lecture_de_TLE.elements_keplerien("données/galileo.txt")['Ascension'][11])

M_galileo = lecture_de_TLE.elements_keplerien("données/galileo.txt")['Anomalie moyenne'][11]

vitesse_angulaire = -7.2921151467e-5

heure_tle = fonctions.tle_to_utc(float(lecture_de_TLE.elements_keplerien("données/galileo.txt")['Heure'][11]))

date_tle = fonctions.t(fonctions.dateJulienne(2021,11,8,heure_tle))


X = []

Y = []

semaine_bis = ['08/11/2021','09/11/2021','10/11/2021','11/11/2021','12/11/2021','13/11/2021','14/11/2021']


semaine = [2459526.5,2459527.5,2459528.5,2459529.5,2459530.5,2459531.5,2459533.49931]

for k in range(1,len(semaine)):
    DATE = np.linspace(fonctions.t(semaine[k-1]),fonctions.t(semaine[k]),num=10000)

    for date in DATE:
        w_galileo = fonctions.var_omega(w0,n0,J2,a_galileo,e_galileo,i_galileo,date,date_tle)

        W_galileo = fonctions.var_OMEGA(W0,n0,J2,a_galileo,e_galileo,i_galileo,date,date_tle)

        M_corr_galileo = fonctions.var_anomalie_moyenne(M_galileo,n0,J2,a_galileo,e_galileo,i_galileo,date,date_tle)

        E_galileo = fonctions.kepler_resolve(M_corr_galileo,e_galileo,1e-6)

        coord_orbital_galileo = fonctions.coord_orbital(a_galileo,e_galileo,E_galileo)

        coord_celeste_galileo = fonctions.coord_celeste(i_galileo,w_galileo,W_galileo,e_galileo,E_galileo,a_galileo)

        coord_cart_galileo = fonctions.coord_cart(i_galileo,w_galileo,W_galileo,e_galileo,E_galileo,a_galileo,date,vitesse_angulaire)

        coord_geodesique_galileo = fonctions.coord_geodesique(i_galileo,w_galileo,W_galileo,E_galileo,date,vitesse_angulaire)

        X.append(np.degrees(coord_geodesique_galileo[0]))
        Y.append(np.degrees(coord_geodesique_galileo[1]))


    world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))

    fig, ax = plt.subplots()
    world.plot(ax=ax)
    plt.scatter(X,Y,s=10,c='red')
    plt.title(f'Trace du satellite Galileo GSAT0209 du {semaine_bis[0]} à {semaine_bis[k]}')
    plt.xlabel('Longitude(deg)')
    plt.ylabel('Latitude(deg)')
    plt.show(block=False)
    plt.pause(3)
    plt.close()


