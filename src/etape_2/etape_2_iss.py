import numpy as np
import fonctions
import fonctions_2
import matplotlib.pyplot as plt
import lecture_de_TLE
import lecture_DORIS


#pense à rajouter la fonction elements_keplerien
# Définition des paramètres du satellite iss BIIF-9

a_iss = fonctions.calcul_a(0.0011264398169714)

i_iss = np.radians(51.6456)

e_iss = 0.0003349

n0 = 0.0011264398169714

J2 = 1.08262652304819194e-03

vitesse_angulaire = -7.2921151467e-5

date = np.linspace(fonctions.t(2459526.5),fonctions.t(2459527.5),num=10000)

for j in range(15):
    azi = []
    ele = []
    w0 = np.radians(lecture_de_TLE.elements_keplerien("données/iss.txt")['Argument du périgée'][0])

    W0 = np.radians(lecture_de_TLE.elements_keplerien("données/iss.txt")['Ascension'][0])


    
    M_iss = M_iss = lecture_de_TLE.elements_keplerien("données/iss.txt")['Anomalie moyenne'][0]

    heure_tle = fonctions.tle_to_utc(float(lecture_de_TLE.elements_keplerien("données/iss.txt")['Heure'][0]))
    date_tle = fonctions.t(fonctions.dateJulienne(2021,11,7,heure_tle))

    name = lecture_DORIS.NAME[j]

    X_sta = lecture_DORIS.X[j]
    Y_sta = lecture_DORIS.Y[j]
    Z_sta = lecture_DORIS.Z[j]

    lam,phi = fonctions_2.coord_geodesique_bis(X_sta,Y_sta,Z_sta)

    for t in date:

        w_iss = fonctions.var_omega(w0,n0,J2,a_iss,e_iss,i_iss,t,date_tle)

        w_iss = fonctions.var_omega(w0,n0,J2,a_iss,e_iss,i_iss,t,date_tle)

        W_iss = fonctions.var_OMEGA(W0,n0,J2,a_iss,e_iss,i_iss,t,date_tle)

        M_corr_iss = fonctions.var_anomalie_moyenne(M_iss,n0,J2,a_iss,e_iss,i_iss,t,date_tle)

        E_iss = fonctions.kepler_resolve(M_corr_iss,e_iss,1e-6)

        coord_cart_iss = fonctions.coord_cart(i_iss,w_iss,W_iss,e_iss,E_iss,a_iss,t,vitesse_angulaire)

        X_sat = coord_cart_iss[0][0]
        Y_sat = coord_cart_iss[1][0]
        Z_sat = coord_cart_iss[2][0]
        if fonctions_2.azimuth_elevation(X_sta,Y_sta,Z_sta,X_sat,Y_sat,Z_sat,lam,phi)[0]>0:
            ele.append(90-np.degrees(fonctions_2.azimuth_elevation(X_sta,Y_sta,Z_sta,X_sat,Y_sat,Z_sat,lam,phi)[0]))
            azi.append(fonctions_2.azimuth_elevation(X_sta,Y_sta,Z_sta,X_sat,Y_sat,Z_sat,lam,phi)[1])
        
    fig = plt.figure()
    ax = fig.add_subplot(projection='polar')
    c = ax.scatter(azi,ele,c='red',s=1) 
    plt.title(f'Skyplot pour la constellation iss depuis la station {name}')
    plt.show(block=False)
    plt.pause(3)
    plt.close()


