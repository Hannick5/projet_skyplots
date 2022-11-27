import numpy as np

def coord_sat_topo(X_sta,Y_sta,Z_sta,X_sat,Y_sat,Z_sat,lam,phi):
    '''Calcul les coordonnées du satellite dans le repère topocentrique

    :param X_sta: X en coordonnées ECEF de la station
    :type X_sta: float
    :param Y_sta: Y en coordonnées ECEF de la station
    :type Y_sta: float
    :param Z_sta: Z en coordonnées ECEF de la station
    :type Z_sta: float
    :param X_sat: X en coordonnées ECEF du satellite
    :type X_sat: float
    :param Y_sat: Y en coordonnées ECEF du satellite
    :type Y_sat: float
    :param Z_sat: Z en coordonnées ECEF du satellite
    :type Z_sat: float
    :param lam: longitude de la station
    :type lam: float
    :param phi: latitude de la station
    :type phi: float

    :return: coordonnées topocentriques du satellite
    :rtype: tuple
    '''

    matrice_de_passage = np.array([[-np.sin(lam),np.cos(lam),0],[-np.cos(lam)*np.sin(phi),-np.sin(lam)*np.sin(phi),np.cos(phi)],[np.cos(lam)*np.cos(phi),np.sin(lam)*np.cos(phi),np.sin(phi)]])
    A = np.array([[X_sat-X_sta],[Y_sat-Y_sta],[Z_sat-Z_sta]])
    coord = np.dot(matrice_de_passage,A)
    X_sat_loc = coord[0][0]
    Y_sat_loc = coord[1][0]
    Z_sat_loc = coord[2][0]
    return (X_sat_loc,Y_sat_loc,Z_sat_loc)


def azimuth_elevation(X_sta,Y_sta,Z_sta,X_sat,Y_sat,Z_sat,lam,phi):
    '''Calcul l'élévation et l'azimuth

    :param X_sta: X en coordonnées ECEF de la station
    :type X_sta: float
    :param Y_sta: Y en coordonnées ECEF de la station
    :type Y_sta: float
    :param Z_sta: Z en coordonnées ECEF de la station
    :type Z_sta: float
    :param X_sat: X en coordonnées ECEF du satellite
    :type X_sat: float
    :param Y_sat: Y en coordonnées ECEF du satellite
    :type Y_sat: float
    :param Z_sat: Z en coordonnées ECEF du satellite
    :type Z_sat: float
    :param lam: longitude de la station
    :type lam: float
    :param phi: latitude de la station
    :type phi: float

    :return: élévation et azimuth
    :rtype: tuple
    '''
    X_sat_loc = coord_sat_topo(X_sta,Y_sta,Z_sta,X_sat,Y_sat,Z_sat,lam,phi)[0]
    Y_sat_loc = coord_sat_topo(X_sta,Y_sta,Z_sta,X_sat,Y_sat,Z_sat,lam,phi)[1]
    Z_sat_loc = coord_sat_topo(X_sta,Y_sta,Z_sta,X_sat,Y_sat,Z_sat,lam,phi)[2]
    
    d = np.sqrt(X_sat_loc**2+Y_sat_loc**2+Z_sat_loc**2)
    d_hori = np.sqrt(X_sat_loc**2+Y_sat_loc**2)
    azimuth = (2*np.arctan(X_sat_loc/(Y_sat_loc +d_hori)))%(2*np.pi)
    elevation = np.arcsin(Z_sat_loc/d)
    return (elevation,azimuth)

def coord_geodesique_bis(X,Y,Z):
    '''Passe les coordonnées cartésiennes en coordonnées géodésiques 

    :param X: X en coordonnées cartésiennes
    :type X: float
    :param Y: Y en coordonnées cartésiennes
    :type Y: float
    :param Z: Z en coordonnées cartésiennes
    :type Z: float

    :return: longitude et latitude
    :rtype: tuple
    '''
    a = 6378137.0
    f = 1/298.257222101
    b = a*(1-f)
    e = np.sqrt((a**2-b**2)/a**2)
    R = np.sqrt(X**2+Y**2+Z**2)
    mu = np.arctan((Z/np.sqrt(X**2+Y**2))*((1-f)+(((e**2)*a)/R)))
    lon = np.arctan2(Y,X)
    lat = np.arctan((Z*(1-f)+(e**2)*a*(np.sin(mu)**3))/((1-f)*(np.sqrt(X**2+Y**2)-(e**2)*a*(np.cos(mu)**3))))
    return (lon,lat)