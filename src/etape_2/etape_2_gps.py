import numpy as np
import fonctions
import fonctions_2
import matplotlib.pyplot as plt
import lecture_de_TLE
import lecture_DORIS

#Définition des paramètres du satellite GPS BIIF-9

a_gps = fonctions.calcul_a(1.45*(10**(-4)))

i_gps = np.radians(53.8532)

e_gps = 0.0063356

n0 = 1.45*(10**(-4))

J2 = 1.08262652304819194e-3

vitesse_angulaire = -7.2921151467e-5

date = np.linspace(fonctions.t(2459526.5),fonctions.t(2459527.5),num=1000)

for j in range(15):
    azi = []
    ele = []
    for i in range(30):
        w0 = np.radians(lecture_de_TLE.elements_keplerien("données/gps-ops.txt")['Argument du périgée'][i])

        W0 = np.radians(lecture_de_TLE.elements_keplerien("données/gps-ops.txt")['Ascension'][i])

        M_gps = lecture_de_TLE.elements_keplerien("données/gps-ops.txt")['Anomalie moyenne'][i]

        heure_tle = fonctions.tle_to_utc(float(lecture_de_TLE.elements_keplerien("données/gps-ops.txt")['Heure'][i]))

        date_tle = fonctions.t(fonctions.dateJulienne(2021,11,8,heure_tle))

        name = lecture_DORIS.NAME[j]

        X_sta = lecture_DORIS.X[j]
        Y_sta = lecture_DORIS.Y[j]
        Z_sta = lecture_DORIS.Z[j]

        lam,phi = fonctions_2.coord_geodesique_bis(X_sta,Y_sta,Z_sta)

        for t in date:

            w_gps = fonctions.var_omega(w0,n0,J2,a_gps,e_gps,i_gps,t,date_tle)

            W_gps = fonctions.var_OMEGA(W0,n0,J2,a_gps,e_gps,i_gps,t,date_tle)

            M_corr_gps = fonctions.var_anomalie_moyenne(M_gps,n0,J2,a_gps,e_gps,i_gps,t,date_tle)

            E_gps = fonctions.kepler_resolve(M_corr_gps,e_gps,1e-6)

            coord_cart_gps = fonctions.coord_cart(i_gps,w_gps,W_gps,e_gps,E_gps,a_gps,t,vitesse_angulaire)

            X_sat = coord_cart_gps[0][0]
            Y_sat = coord_cart_gps[1][0]
            Z_sat = coord_cart_gps[2][0]
            if fonctions_2.azimuth_elevation(X_sta,Y_sta,Z_sta,X_sat,Y_sat,Z_sat,lam,phi)[0]>0:
                ele.append(90-np.degrees(fonctions_2.azimuth_elevation(X_sta,Y_sta,Z_sta,X_sat,Y_sat,Z_sat,lam,phi)[0]))
                azi.append(fonctions_2.azimuth_elevation(X_sta,Y_sta,Z_sta,X_sat,Y_sat,Z_sat,lam,phi)[1])
            
    fig = plt.figure()
    ax = fig.add_subplot(projection='polar')
    ax.set_theta_zero_location("N")
    c = ax.scatter(azi,ele,c='red',s=1) 
    plt.title(f'Skyplot pour la constellation GPS depuis la station {name}')
    plt.show(block=False)
    plt.pause(3)
    plt.close()


