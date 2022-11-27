from math import radians
import fonctions
import lecture_de_TLE
import numpy as np

def main():
    """Cette fonction a pour but de tester les fonctions qu'on a créees à l'étape 1"""

    # Calcul de a

    print("La valeur de a est: "+str(fonctions.calcul_a(1.45*(10**(-4))))+" mètres\n")

    # Pour tester les fonctions suivantes on va considérer le satellite GPS BIIR-2 (PRN 13) dans le fichier gps-ops.txt

    satellite = lecture_de_TLE.elements_keplerien("données/gps-ops.txt")
    
    # On teste la fonction tle to utc et on vérifie par ailleurs que l'appel d'éléments dans le TLE fonctionne
    
    print("L'heure utc est " +str(fonctions.tle_to_utc(float(satellite['Heure'][0])))+"pour le satellite GPS BIIR-2.\n")

    # On a bien le même résultat donc l'appel d'éléments fonctionne 

    # Essayons avec l'exemple 8.6 du poly qui doit donner 6:00:15.79

    print("Expected : "+str((6,0,15.79))) 
    print("Get : "+str(fonctions.tle_to_utc(0.25018279))+"\n")

    # On teste la résolution de l'équation de Kepler
    M = np.radians(274.5785)
    e = 0.0002250
    a = 6971515
    print("La résolution de l'équation de Kepler pour M = "+str(M)+" donne E = "+str(fonctions.kepler_resolve(M,e,1e-6))+"\n")

    # On n'a pas la même valeur que dans le poly mais l'itération fonctionne et nous donne une valeur de E proche de M puisque e est petit

    # On teste la fonction qui calcule le rayon orbital et l'anomalie vraie

    print("le rayon orbital est de : "+str(fonctions.rayon_orbital(a,e,np.radians(274.566)))+" mètres\n")
    r = fonctions.rayon_orbital(a,e,np.radians(274.566))
    E = np.radians(274.566)
    print("Expected : v = "+str(1.49)) #calculatrice
    print("Get : v = "+str(fonctions.anomalie_vraie(e,a,E,r))+"\n")
    v = fonctions.anomalie_vraie(e,a,E,r)

    # On teste la fonction qui calcule les coordonnées orbitales

    print("Les coordonnées orbitales sont : "+str(fonctions.coord_orbital(a,e,E))+"\n")

    # On reprend le cas du satellite 
    e_sat = 0.0054755
    M_sat = satellite['Anomalie moyenne'][0]
    E_sat = fonctions.kepler_resolve(M_sat,e_sat,1e-6)
    i_sat = np.radians(55.4827)
    a_sat = fonctions.calcul_a(1.45*(10**(-4)))
    r_sat = fonctions.rayon_orbital(a_sat,e_sat,E_sat)
    v_sat = fonctions.anomalie_vraie(e_sat,a_sat,E_sat,r_sat)
    
    print("Les coordonnées orbitales sont : "+str(fonctions.coord_orbital(a_sat,e_sat,E_sat))+"\n")
    
    # Ceci nous renvoie bien un tuple contenant (x,y)

    #Test du passage en coordonnées céleste
    w_sat = np.radians(satellite['Argument du périgée'][0])
    W_sat = np.radians(satellite['Ascension'][0])
    print('Les coordonnées céleste sont : \n'+str(fonctions.coord_celeste(i_sat,w_sat,W_sat,e_sat,E_sat,a_sat))+"\n")

    # On obtient un ndarray contenant les coordonnées célestes

    #Test du passage en coordonnées cartesienne
    date = fonctions.t(fonctions.dateJulienne(2021,11,8,0)) #On définit une date fixe à laquelle appliqué la correction ici 08/11/2021 00h00
    vitesse_angulaire = -7.2921151467e-5
    print("Les coordonnées cartésiennes sont : \n"+str(fonctions.coord_cart(i_sat,w_sat,W_sat,e_sat,E_sat,a_sat,date,vitesse_angulaire))+"\n")

    # On obtient un ndarray contenant les coordonnées cartésiennes

    # Test passage en coordonnées géodésiques

    print("Les coordonnées géodésiques sont : \n" +str(fonctions.coord_geodesique(i_sat,w_sat,W_sat,E_sat,date,vitesse_angulaire)))

    # On obtient un tuple avec la longitude et la latitude en radian 

    




    
if __name__ == "__main__":
    main()
