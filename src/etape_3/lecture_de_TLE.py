from pathlib import Path
import pandas as pd
import math


# Lecture des données identique à celle de lecture_des_tle.py
working_directory = Path(__file__).absolute().parent


def prepend(list, str): 
    '''Insérer une chaîne de caractère au début de tous les éléments d'une liste

    :param list: une liste
    :type M: list
    :param str: une chaîne de caractère
    :type e: string

    :return: liste des valeurs modifiées
    :rtype: list 
    '''
     
    str += '{0}'
    list = [str.format(i) for i in list]
    return(list)


def elements_keplerien(fichier):
    '''Créer un tableau pandas rempli d'éléments képlériens format [ascension droite argument du périgée anormalie moyenne]

    :rtype: pandas.core.frame.DataFrame
    :return: tableau pandas d'éléments képlériens
    '''

    omega=[]
    petit_omega =[]
    M= []
    nom_sat_beta = []
    date =[]
   
    with open(working_directory/fichier,"r") as f:
        # Une liste contenant chaque ligne du fichier comme un élément de liste.
        lines = f.readlines()

        # Parcourir les lignes du fichier TLE
        for l in lines:
            # Parcourir les valeurs des lignes d'indice 2 du fichier TLE
            if l[0] == '2': 
                omega.append(l.split()[3])
                petit_omega.append(l.split()[5])
                M.append(l.split()[6])

            # Parcourir les valeurs des lignes d'indice 1 du fichier TLE
            if l[0] == '1':
                date.append(l.split()[3])

            # Obtenir les noms des satellites; ils commencent soit par C, soit par G
            if l[0] == 'G' or l[0] == 'C' or l[0] == 'I':
                nom_sat_beta.append(l)

    # Transformer les valeurs du type string en float
    omega = [ float(x) for x in omega ]
    petit_omega = [ float(x) for x in petit_omega ]
    M = [ float(x) for x in M ]  
   
   
    annee_jour = []
    heure = []

    # Séparer les valeurs avant et après la virgule des éléments 'date'
    for i in range(len(date)):
        annee_jour.append(date[i].split('.',1)[0])
        heure.append(date[i].split('.',1)[1])
   
    # Ajouter 0. à chaque élément de la liste heure
    heure = prepend(heure, '0.')
   
    # Obtenir les deux derniers chiffres avant la virgule qui représentent le jour
    jour = [e[2:] for e in annee_jour]
   
    # Transformer les données de M de degré en rad  
    for i in range(len(M)):
        M[i] =(M[i]*math.pi)/180
       
    # Enlever les /n dans la liste nom_sat
    nom_sat = [s.replace("\n", "") for s in nom_sat_beta]
   
    # Créer un tableau avec pandas
    mydataset = {
  'Ascension': omega, 'Argument du périgée' : petit_omega,
  'Anomalie moyenne' : M, 'Heure' : heure, 'Jour' : jour
}
    myvar = pd.DataFrame(mydataset, index = nom_sat)
   
    # Afficher toutes les colonnes du tableau, car certaines sont cachées
    pd.set_option('display.max_columns', None)
    return myvar