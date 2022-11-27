import numpy as np
import fonctions
import fonctions_2
import lecture_de_TLE
import matplotlib.pyplot as plt

# Paramètres ISS

a_iss = fonctions.calcul_a(0.0011264398169714)

i_iss = np.radians(51.6456)

e_iss = 0.0003349

n0_iss = 0.0011264398169714

J2 = 1.08262652304819194e-03

vitesse_angulaire = -7.2921151467e-5

date = np.linspace(fonctions.t(2459526.5),fonctions.t(2459527.5),num=10000)

w0_iss = np.radians(lecture_de_TLE.elements_keplerien("données/iss.txt")['Argument du périgée'][0])

W0_iss = np.radians(lecture_de_TLE.elements_keplerien("données/iss.txt")['Ascension'][0])



M_iss = M_iss = lecture_de_TLE.elements_keplerien("données/iss.txt")['Anomalie moyenne'][0]

heure_tle_iss = fonctions.tle_to_utc(float(lecture_de_TLE.elements_keplerien("données/iss.txt")['Heure'][0]))
date_tle_iss = fonctions.t(fonctions.dateJulienne(2021,11,7,heure_tle_iss))

# Paramètres GPS

a_gps = fonctions.calcul_a(1.45*(10**(-4)))

i_gps = np.radians(53.8532)

e_gps = 0.0063356

n0 = 1.45*(10**(-4))

NAME = lecture_de_TLE.elements_keplerien("données/gps-ops.txt").index.tolist()


for i in range(30):
    azi = []
    ele = []

    w0 = np.radians(lecture_de_TLE.elements_keplerien("données/gps-ops.txt")['Argument du périgée'][i])

    W0 = np.radians(lecture_de_TLE.elements_keplerien("données/gps-ops.txt")['Ascension'][i])

    M_gps = lecture_de_TLE.elements_keplerien("données/gps-ops.txt")['Anomalie moyenne'][i]

    heure_tle = fonctions.tle_to_utc(float(lecture_de_TLE.elements_keplerien("données/gps-ops.txt")['Heure'][i]))

    date_tle = fonctions.t(fonctions.dateJulienne(2021,11,8,heure_tle))
        
    for t in date:

        w_gps = fonctions.var_omega(w0,n0,J2,a_gps,e_gps,i_gps,t,date_tle)

        W_gps = fonctions.var_OMEGA(W0,n0,J2,a_gps,e_gps,i_gps,t,date_tle)

        M_corr_gps = fonctions.var_anomalie_moyenne(M_gps,n0,J2,a_gps,e_gps,i_gps,t,date_tle)

        E_gps = fonctions.kepler_resolve(M_corr_gps,e_gps,1e-6)

        coord_cart_gps = fonctions.coord_cart(i_gps,w_gps,W_gps,e_gps,E_gps,a_gps,t,vitesse_angulaire)

        X_sat = coord_cart_gps[0][0]
        Y_sat = coord_cart_gps[1][0]
        Z_sat = coord_cart_gps[2][0]


        w_iss = fonctions.var_omega(w0_iss,n0_iss,J2,a_iss,e_iss,i_iss,t,date_tle_iss)

        w_iss = fonctions.var_omega(w0_iss,n0_iss,J2,a_iss,e_iss,i_iss,t,date_tle_iss)

        W_iss = fonctions.var_OMEGA(W0_iss,n0_iss,J2,a_iss,e_iss,i_iss,t,date_tle_iss)

        M_corr_iss = fonctions.var_anomalie_moyenne(M_iss,n0_iss,J2,a_iss,e_iss,i_iss,t,date_tle_iss)

        E_iss = fonctions.kepler_resolve(M_corr_iss,e_iss,1e-6)

        coord_cart_iss = fonctions.coord_cart(i_iss,w_iss,W_iss,e_iss,E_iss,a_iss,t,vitesse_angulaire)

        X_antenne = coord_cart_iss[0][0]
        Y_antenne = coord_cart_iss[1][0]
        Z_antenne = coord_cart_iss[2][0]

        lam,phi = fonctions_2.coord_geodesique_bis(X_antenne,Y_antenne,Z_antenne)

        if fonctions_2.azimuth_elevation(X_antenne,Y_antenne,Z_antenne,X_sat,Y_sat,Z_sat,lam,phi)[0]>0:
                ele.append(90-np.degrees(fonctions_2.azimuth_elevation(X_antenne,Y_antenne,Z_antenne,X_sat,Y_sat,Z_sat,lam,phi)[0]))
                azi.append(fonctions_2.azimuth_elevation(X_antenne,Y_antenne,Z_antenne,X_sat,Y_sat,Z_sat,lam,phi)[1])

    fig = plt.figure()
    ax = fig.add_subplot(projection='polar')
    ax.set_theta_zero_location("N")
    c = ax.scatter(azi,ele,c='red',s=1) 
    plt.title(f"Skyplot pour {NAME[i].replace(' ','')} depuis l'antenne sur l'ISS")
    plt.show(block=False)
    plt.pause(3)
    plt.close()



