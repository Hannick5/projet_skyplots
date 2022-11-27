import numpy as np
import fonctions
import fonctions_2
import matplotlib.pyplot as plt
import lecture_de_TLE
import lecture_DORIS

# Définition des paramètres du satellite glo BIIF-9

a_glo = fonctions.calcul_a(1.45*(10**(-4)))

i_glo = np.radians(65.9702)

e_glo = 0.0025200

J2 = 1.08262652304819194e-3

n0 = 1.45*(10**(-4))

vitesse_angulaire = -7.2921151467e-5

date = np.linspace(fonctions.t(2459526.5),fonctions.t(2459527.5),num=1000)

azi = []
ele = []

for j in range(15):
    azi = []
    ele = []
    for i in range(27):
        w0 = np.radians(lecture_de_TLE.elements_keplerien("données/glo-ops.txt")['Argument du périgée'][i])

        W0 = np.radians(lecture_de_TLE.elements_keplerien("données/glo-ops.txt")['Ascension'][i])

        M_glo = lecture_de_TLE.elements_keplerien("données/glo-ops.txt")['Anomalie moyenne'][i]

        heure_tle = fonctions.tle_to_utc(float(lecture_de_TLE.elements_keplerien("données/glo-ops.txt")['Heure'][i]))

        date_tle = fonctions.t(fonctions.dateJulienne(2021,11,8,heure_tle))

        name = lecture_DORIS.NAME[j]



        X_sta = lecture_DORIS.X[j]
        Y_sta = lecture_DORIS.Y[j]
        Z_sta = lecture_DORIS.Z[j]

        lam,phi = fonctions_2.coord_geodesique_bis(X_sta,Y_sta,Z_sta)

        for t in date:
            w_glo = fonctions.var_omega(w0,n0,J2,a_glo,e_glo,i_glo,t,date_tle)

            w_glo = fonctions.var_omega(w0,n0,J2,a_glo,e_glo,i_glo,t,date_tle)

            W_glo = fonctions.var_OMEGA(W0,n0,J2,a_glo,e_glo,i_glo,t,date_tle)

            M_corr_glo = fonctions.var_anomalie_moyenne(M_glo,n0,J2,a_glo,e_glo,i_glo,t,date_tle)

            E_glo = fonctions.kepler_resolve(M_corr_glo,e_glo,1e-6)

            coord_cart_glo = fonctions.coord_cart(i_glo,w_glo,W_glo,e_glo,E_glo,a_glo,t,vitesse_angulaire)

            X_sat = coord_cart_glo[0][0]
            Y_sat = coord_cart_glo[1][0]
            Z_sat = coord_cart_glo[2][0]
            if fonctions_2.azimuth_elevation(X_sta,Y_sta,Z_sta,X_sat,Y_sat,Z_sat,lam,phi)[0]>0:
                ele.append(90-np.degrees(fonctions_2.azimuth_elevation(X_sta,Y_sta,Z_sta,X_sat,Y_sat,Z_sat,lam,phi)[0]))
                azi.append(fonctions_2.azimuth_elevation(X_sta,Y_sta,Z_sta,X_sat,Y_sat,Z_sat,lam,phi)[1])
            
    fig = plt.figure()
    ax = fig.add_subplot(projection='polar')
    c = ax.scatter(azi,ele,c='red',s=1) 
    plt.title(f'Skyplot pour la constellation glo depuis la station {name}')
    plt.show(block=False)
    plt.pause(3)
    plt.close()


