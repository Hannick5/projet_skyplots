# Lecture du fichier DORIS

from pathlib import Path

working_directory = Path(__file__).absolute().parent

X=[]
Y=[]
Z=[]
NAME = []

with open(working_directory/'données/DORIS_stations.txt','r') as f:

    # Parcourir les lignes du fichier en sautant la toute première ligne
    lines = f.readlines()[1:]

    for l in lines:
        NAME.append(l.split()[0])
        X.append(l.split()[1])
        Y.append(l.split()[2])
        Z.append(l.split()[3])

# Transformer les valeurs de X Y Z en float
X = [float(x) for x in X]
Y = [float(x) for x in Y]
Z = [float(x) for x in Z]

