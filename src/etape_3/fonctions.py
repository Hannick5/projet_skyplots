import pandas as pd 
import numpy as np


# Calcul du demi-grand axe de l'orbite du satellite a
def calcul_a(n):
    '''Calcule le paramètre a en fonction du moyen mouvement

    :param n: [moyen mouvement ramené en rad/s]
    :type n: [float]

    :return: le demi-grand de l'orbite du satellite
    :rtype: float
    '''

    GM = 3.986005*(10**14)
    a = (GM/(n**2))**(1/3)
    return a
    
# On trouve a = 26664 km

# On définit t
def dateJulienne(annee,mois,jour,heure):
    '''Calcule la date julienne à partir de la date format annee mois jour heure

    :param annee: année
    :type annee: int
    :param mois: mois
    :type mois: int
    :param jour: jour
    :type jour: int
    :param heure: heure utc
    :type heure: int

    :return: date julienne
    :rtype: int
    '''

    if mois<=2:
        a=annee-1
        m=mois+12
    else:
        a=annee
        m=mois
    b=int(a/400)-int(a/100)
    DJ = int(365.25*a)+int(30.6001*(m+1))+b+1720996.5+jour+heure*(1/24)
    return DJ



def t(DJ):
    '''Calcul date julienne modifiée avec la date julienne

    :param DJ: date julienne
    :type DJ: int

    :return: date julienne modifiée
    :rtype: float
    '''
    
    return (DJ-2400000.5)*86400

# On code une fonction qui va convertir les heures du TLE en heure UTC
def tle_to_utc(heure_tle):
    heure = heure_tle*24
    return np.floor(heure)

# Calcul de la variation d’anomalie moyenne depuis l’instant de référence du TLE
def var_anomalie_moyenne(M,n0,J2,a,e,i,t,t_tle):
    
    R = 6371000
    deltan = (3/4)*n0*J2*((R/a)**2)*(1/((1-e*e)**(3/2)))*(3*(np.cos(i)**2)-1)
    delta_M = (n0+deltan)*(t-t_tle)
    M = M + delta_M
    return M

# Calcul de la variation pour l'Ascension
def var_OMEGA(W0,n0,J2,a,e,i,t,t0):
    R = 6371000
    W_dot = (-3/2)*n0*J2*((R/a)**2)*(1/((1-e*e)**2))*np.cos(i)
    W = W0 + W_dot*(t-t0)
    return W

# Calcul de la variation pour l'argument

def var_omega(w0,n0,J2,a,e,i,t,t0):
    R = 6371000
    w_dot = (3/4)*n0*J2*((R/a)**2)*(1/((1-e*e)**2))*(5*(np.cos(i)**2-1))
    w = w0 + w_dot*(t-t0)
    return w


# Résolution équation de Kepler (passage de l'anomalie moyenne à l'anomalie excentrique)

def kepler_resolve(M,e,epsilon):
    '''Calcul la valeur de E avec la méthode de résolution de l'équation de Kepler

    :param M: anomalie moyenne du satellite en radian
    :type M: list
    :param e: excentricité du satellite
    :type e: float
    :param epsilon: précision/erreur
    :type epsilon: float

    :return: la valeur de l'anomalie excentrique
    :rtype: float
    '''

    E = M
    seuil = 1e6
    while seuil>epsilon:
        nouveau_E = M + e*np.sin(E)
        seuil = abs(nouveau_E-E)
        E = nouveau_E
    return E

# Calcul du rayon orbital
def rayon_orbital(a,e,E):
    '''Calcul le rayon orbital pour un satellite donné

    :param a: demi-grand axe
    :type a: float
    :param e: excentricité
    :type e: float
    :param E: anomalie excentrique
    :type E: float

    :return: le rayon orbital
    :rtype: float
    '''
    r = a*(1-e*np.cos(E))
    return r


# Calcul de l'anomalie vraie
def anomalie_vraie(e,a,E,r):
    '''Calcul l'anomalie vraie

    :param a: demi grand axe
    :type a: float
    :param e: excentricité
    :type e: float
    :param E: anomalie excentrique
    :type E: float
    :param r: rayon orbital
    :type r: float

    :return: l'anomalie vraie
    :rtype: float
    '''
    v = np.arccos((a/r)*(np.cos(E)-e)) #E est en radian comme on a convertis M en radian
    return v



# Coordonnées dans le plan orbital
def coord_orbital(a,e,E):
    '''Calcul les coordonnées dans le plan orbital

    :param e: excentricité
    :type e: float
    :param E: anomalie excentrique
    :type E: float

    :return: couple de coordonnées
    :rtype: tuple
    '''

    x = a*(np.cos(E)-e)
    y = a*np.sqrt(1-e*e)*np.sin(E)
    return (x,y)


# Coordonnées dans le plan céleste
def coord_celeste(i,w,W,e,E,a):
    '''Passage des coordonnées orbitales en coordonnée célestes

    :param i: inclinaison
    :type i: float
    :param w: argument du périgée
    :type w: float
    :param W: ascension droite
    :type W: float
    :param r: rayon orbital
    :type r: float
    :param v: anomalie vraie
    :type v: float

    :return: coordonnées célestes
    :rtype: tuple
    '''

    R_inclinaison = np.array([[1,0,0],[0,np.cos(i),-np.sin(i)],[0,np.sin(i),np.cos(i)]])
    R_argument = np.array([[np.cos(w),-np.sin(w),0],[np.sin(w),np.cos(w),0],[0,0,1]])
    R_ascension = np.array([[np.cos(W),-np.sin(W),0],[np.sin(W),np.cos(W),0],[0,0,1]])
    (x,y) = coord_orbital(a,e,E)
    coord_orb = np.array([[x],[y],[0]])
    coord = R_ascension.dot(R_inclinaison).dot(R_argument.dot(coord_orb))

    return coord

# Coordonnées cartésiennes
def coord_cart(i,w,W,e,E,a,t,vitesse_angulaire):
    '''Passage des coordonnées célestes aux coordonnées cartésiennes ECEF

    :param i: inclinaison
    :type i: float
    :param w: argument du périgée
    :type w: float
    :param W: ascension droite
    :type W: float
    :param r: rayon orbital
    :type r: float
    :param v: anomalie vraie
    :type v: float
    :param t: date
    :type t: int
    :param vitesse_angulaire: vitesse angulaire de rotation de la Terre
    :type vitesse_angulaire: float

    :return: coordonnées cartésiennes
    :rtype: tuple
    '''
    
    matrice_passage = np.array([[np.cos(vitesse_angulaire*t),-np.sin(vitesse_angulaire*t),0],[np.sin(vitesse_angulaire*t),np.cos(vitesse_angulaire*t),0],[0,0,1]])
    coord_car = np.dot(matrice_passage,coord_celeste(i,w,W,e,E,a))
    return coord_car

# Coordonnées géodésiques ECEF
def coord_geodesique(i,w,W,E,t,vitesse_angulaire):
    '''Passage des coordonnées cartésiennes ECEF aux coordonnées géodésiques

    :param i: inclinaison
    :type i: float
    :param w: argument du périgée
    :type w: float
    :param W: ascension droite
    :type W: float
    :param r: rayon orbital
    :type r: float
    :param vitesse_angulaire: vitesse angulaire de rotation de la Terre
    :type vitesse_angulaire: float
  
    :return: coordonnées géodésiques
    :rtype: tuple
    '''

    a = 6378137.0
    f = 1/298.257222101
    b = a*(1-f)
    e = np.sqrt((a**2-b**2)/a**2)
    coord_car = coord_cart(i,w,W,e,E,a,t,vitesse_angulaire)
    X = coord_car[0][0]
    Y = coord_car[1][0]
    Z = coord_car[2][0]
    R = np.sqrt(X**2+Y**2+Z**2)
    mu = np.arctan((Z/np.sqrt(X**2+Y**2))*((1-f)+(((e**2)*a)/R)))
    lon = np.arctan2(Y,X)
    lat = np.arctan((Z*(1-f)+(e**2)*a*(np.sin(mu)**3))/((1-f)*(np.sqrt(X**2+Y**2)-(e**2)*a*(np.cos(mu)**3))))
    return (lon,lat)








