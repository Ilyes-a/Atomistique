import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle


def atomistique(X, action):
    """Les modules matplotlib (plt) et numpy (np) doivent etre importes
    Affiche selon la demande :
    - la configuration electronique complete
    - la configuration en notation allégée, exprimée en fonction du gaz noble le plus proche
    - la configuration de valence
    - la représentation en cases quantiques

    Parametres:
        X (int/str) : entite chimique
        action (str) : configuration demandee qui peut prendre les valeurs :
                           - "structureComplete"
                           - "structureAllegee"
                           - "configurationValence"
                           - "casesQuantiques"
    """

    # Ordre de remplissage conforme  à la règle de Klechowski
    ordre = '1s2s2p3s3p4s3d4p5s4d5p6s4f5d6p7s5f6d7p'

    # Dictionnaire contenant le nombre maximum d'électrons par sous-couche
    mSc = {'s': 2, 'p': 6, 'd': 10, 'f': 14}

    # Tuple des gaz nobles
    gazNobles = ('He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn')

    # Tuple contenant tous les éléments chimiques
    elements = ('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',
                'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr',
                'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br',
                'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',
                'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La',
                'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',
                'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',
                'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
                'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md',
                'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',
                'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og')

    # Exceptions à la règle de Klechkowski
    exceptions = ('Cr', 'Cu', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', 'Ag', 'Pt', 'Au',
                  'La', 'Ce', 'Gd', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Cm')

    # Définition des fonctions

    def nombreElectrons(X):
        """Calcule le nombre d'electrons

        Parametres:
            X (int) : entite chimique
            X (str) : entite chimique

        Retours:
            tuple(symbole de l'entite, nombre d'electrons, element)
        """

        ################################################### X est un entier ###################################################

        if type(X) is int:
            if X > 0 and X < 119:  # nombre d'electrons compris entre 1 et 118
                # renvoie l'element et le numero atomique de l'element/nombre d'electrons
                return elements[X-1], X, elements[X-1]
            else:
                # Erreur : X n'est pas dans la plage d'entiers definie
                return "Error : 'X' out of range"

        ############################################ X est une chaine de caracteres ############################################

        elif type(X) is str:

            # conversion de la chaine de caracteres en liste des caracteres
            tab_X = list(str.strip(X))

            # definition des variables
            pos_elt1 = 1
            pos_elt2 = 1
            pos_elt3 = 1
            pos_elt4 = 1
            pos_elt5 = 1
            pos_digit = 0
            pos_plusOuMoins = 0

            # cas ou X est un ion -------------------------------------------------------------------------------
            if any(itm == "+" or itm == "-" for itm in tab_X):

                # ----------- cas ou X est note avec un chiffre et un signe --------------
                if any(n.isdigit() for n in X):  # verifie s'il y a un chiffre dans X
                    for compteur1 in tab_X:  # parcours de toute la liste
                        if compteur1.isdigit():  # recherche du chiffre dans la liste
                            # separation de la liste en deux parties
                            tab_elt1 = tab_X[:pos_digit]  # partie element
                            tab_charge1 = tab_X[pos_digit:]  # partie charge
                            if len(tab_elt1) == 1:
                                elt1 = tab_elt1[0]
                            elif len(tab_elt1) == 2:
                                elt1 = tab_elt1[0]+tab_elt1[1]
                            for compteur2 in elements:  # parcours de tous les elements du tuple
                                if compteur2 == elt1:  # l'element correspond a elt1
                                    if tab_charge1[1] == "+":  # cas cation
                                        # renvoie le symbole de l'entite, le nombre d'electrons et l'element
                                        return X, pos_elt1 - int(tab_charge1[0]), elt1
                                    elif tab_charge1[1] == "-":
                                        # renvoie le symbole de l'entite, le nombre d'electrons et l'element
                                        return X, pos_elt1 + int(tab_charge1[0]), elt1
                                else:
                                    pos_elt1 += 1
                        else:
                            pos_digit += 1

                # --------------- cas ou X est note avec un seul signe -------------------
                elif tab_X.count("+") == 1 or tab_X.count("-") == 1:
                    for compteur3 in elements:  # parcours de tous les elements de la liste
                        if compteur3 == tab_X[0]:  # element a une seule lettre
                            if tab_X[1] == "+":  # cas cation
                                # renvoie le symbole de l'entite, le nombre d'electrons et l'element
                                return X, pos_elt2 - 1, tab_X[0]
                            elif tab_X[1] == "-":  # cas anion
                                # renvoie le symbole de l'entite, le nombre d'electrons et l'element
                                return X, pos_elt2 + 1, tab_X[0]
                        else:
                            pos_elt2 += 1

                        # element a deux  lettres
                        if compteur3 == tab_X[0]+tab_X[1]:
                            if tab_X[2] == "+":  # cas cation
                                # renvoie le symbole de l'entite, le nombre d'electrons et l'element
                                return X, pos_elt3 - 1, tab_X[0]+tab_X[1]
                            elif tab_X[2] == "-":  # cas anion
                                # renvoie le symbole de l'entite, le nombre d'electrons et l'element
                                return X, pos_elt3 + 1, tab_X[0]+tab_X[1]
                        else:
                            pos_elt3 += 1

                # -------------- cas ou X est note avec plusieurs signes -----------------
                else:
                    for x in tab_X:  # parcours de toute la liste
                        if x == "+" or x == "-":  # recherche du premier signe dans la liste
                            # separation de la liste en deux parties
                            # partie element
                            tab_elt2 = tab_X[:pos_plusOuMoins]
                            # partie charge
                            tab_charge2 = tab_X[pos_plusOuMoins:]
                            if len(tab_elt2) == 1:
                                elt2 = tab_elt2[0]
                            elif len(tab_elt2) == 2:
                                elt2 = tab_elt2[0]+tab_elt2[1]
                            for j in elements:  # parcours de tous les elements du tuple
                                if j == elt2:  # l'element correspond a elt2
                                    if tab_charge2[1] == "+":  # cas cation
                                        # renvoie le symbole de l'entite, le nombre d'electrons et l'element
                                        return X, pos_elt4 - len(tab_charge2), elt2
                                    elif tab_charge2[1] == "-":  # cas anion
                                        # renvoie le symbole de l'entite, le nombre d'electrons et l'element
                                        return X, pos_elt4 + len(tab_charge2), elt2
                                else:
                                    pos_elt4 += 1
                        else:
                            pos_plusOuMoins += 1

            # cas ou X est un atome -----------------------------------------------------------------------------
            else:
                for compteur4 in elements:  # parcours de tous les elements du tuple
                    if compteur4 == X:  # l'element correspond a X
                        # renvoie l'element et le numero atomique de l'element/nombre d'electrons
                        return X, pos_elt5, X
                    else:
                        pos_elt5 += 1

    def structureComplete(X):
        """Affiche la configuration electronique complete

        Parametres:
            X (int) : entite chimique
            X (str) : entite chimique

        Retours:
            sc_excep (str) : structure complete faisant exception a la regle de Klechkowski
            sc_verif (str) : structure complete verifiant la regle de Klechkowski
        """

        # definition des variables
        division = list(str.strip(ordre))
        ordreBase = []
        structure_lst = []
        structure_lst2 = []
        structure_str = ""
        x_restant = nombreElectrons(X)[1]
        compteur = 0

        # separation de la chaine de caracteres "ordre" en liste de listes
        for k in range(0, len(division), 2):
            ordreBase.append([division[k], division[k+1]])

        # ajout du nombre d'electrons maximum par sous-couches
        for e in ordreBase:
            e[0] = int(e[0])
            if e[1] == "s":
                e.append(mSc["s"])
            elif e[1] == "p":
                e.append(mSc["p"])
            elif e[1] == "d":
                e.append(mSc["d"])
            elif e[1] == "f":
                e.append(mSc["f"])

        # calcul et affichage de la structure
        while x_restant > 0:
            n = ordreBase[compteur][0]
            nom = ordreBase[compteur][1]
            ne = ordreBase[compteur][2]

            if x_restant < ne:
                nmin = x_restant
            else:
                nmin = ne

            structure_lst.append([n, nom, nmin])
            structure_lst2.append(n)
            structure_lst2.append(nom)
            structure_lst2.append(nmin)
            structure_str = structure_str + str(n) + nom + str(nmin) + " "
            x_restant = x_restant - nmin
            compteur += 1

        # renvoi de la structure complete en fonction du respect de la regle de Klechkowski

        # excetion a la regle de Klechkowski
        if nombreElectrons(X)[2] in exceptions:
            if type(X) is int:
                return elements[X-1] + "*" + " : " + structure_str, structure_lst, structure_lst2
            elif type(X) is str:
                return X + "*" + " : " + structure_str, structure_lst, structure_str, structure_lst2
        # verifie la regle de Klechkowski
        else:
            if type(X) is int:
                return elements[X-1] + " : " + structure_str, structure_lst, structure_lst2
            elif type(X) is str:
                return X + " : " + structure_str, structure_lst, structure_lst2

    def structureAllegee(X):
        """Affiche la configuration en notation allegee

        Parametres:
            X (int) : entite chimique
            X (str) : entite chimique

        Retours:
            sc_finale_str (str) : structure allegee en fonction du gaz noble le plus proche
            structureComplete(X)[0] (str) : structure allegee (= structure complete)
        """

        # definition des variables
        sc_X = structureComplete(X)[1]
        sc_gN = []
        sc_inter = []
        sc_finale = []
        sc_finale_str = ""

        # construction des listes intermediaires
        for k in gazNobles:
            sc_gN.append(structureComplete(k)[1])

        sc_gN_inverse = sc_gN[::-1]
        gazNobles_inverse = gazNobles[::-1]

        # fonction principale
        if len(sc_X) >= 2:  # Li (Z = 3) --> Og (Z = 118)
            while len(sc_X) > 0:
                for i in range(6):
                    for index, j in enumerate(sc_gN_inverse):
                        if j == sc_X:
                            sc_finale.append([gazNobles_inverse[index]])
                            sc_finale.extend(sc_inter)
                            if type(X) is int:
                                sc_finale_str = sc_finale_str + \
                                    nombreElectrons(X)[2] + " :"
                            else:
                                sc_finale_str = sc_finale_str + X + " :"
                            for p in range(0, 1, 3):
                                for n in sc_finale:
                                    sc_finale_str = sc_finale_str + " "
                                    for m in n:
                                        sc_finale_str = sc_finale_str + str(m)
                            return sc_finale_str
                p = sc_X.pop()
                sc_inter.insert(0, p)
        else:  # H (Z = 1) + He (Z = 2)
            return structureComplete(X)[0]

    def configurationValence(X):
        """Affiche la configuration de Valence

        Parametres:
            X (int) : entite chimique
            X (str) : entite chimique

        Retours:
            valence_str (str) : electrons de nombre quantique principal le plus élevé
                                et ceux appartenants aux sous-couches partiellement remplies
        """

        # definition des variables
        maxi = 0
        maxi_liste = []
        part_remplie = []
        key_X = structureComplete(X)[1][-1][1]
        value_X = structureComplete(X)[1][-1][2]
        valence_str = ""

        # recherche du nombre quantique principal le plus eleve
        for k in structureComplete(X)[1]:
            if k[0] > maxi:
                maxi = k[0]
                maxi_liste = [k]
            elif k[0] == maxi:
                maxi_liste.append(k)

        # retour de fonction
        if type(X) is int:
            valence_str = valence_str + "Electrons de valence de " + \
                nombreElectrons(X)[2] + " :"
        else:
            valence_str = valence_str + "Electrons de valence de " + X + " :"

        for i in maxi_liste:
            valence_str = valence_str + " "
            for j in i:
                valence_str = valence_str + str(j)

        valence_str = valence_str + " "

        for n in part_remplie:
            if len(structureComplete(X)[1]) > 1:
                valence_str = valence_str + str(n)
        return valence_str, maxi

    def casesQuantiques(X):
        """Affiche la representation en cases quantiques

        Parametres:
            X (int) : entite chimique
            X (str) : entite chimique

        Retours:
            figure matplotlib
        """

        # definition des variables
        fig = plt.figure(figsize=(12, 7))
        ax = fig.add_subplot()
        cq_max = configurationValence(X)[1]
        l = 25
        h = 2 * cq_max + 1
        h_quantique = 2 * cq_max - 1
        h_case_quantique = 2 * cq_max
        cc = structureComplete(X)[2]
        title = "Structure de " + structureAllegee(X)
        sc = "Structure complète de " + structureComplete(X)[0]

        # construction du rectangle principal
        rect = Rectangle((0, 0), l, h, fill=False, color="black")
        ax.add_patch(rect)

        # construction de la legende
        plt.text(2, 2 * cq_max + 1.5, title, fontsize=10,
                 weight="bold", color="orange")
        plt.text(2, -1, sc, fontsize=10, weight="bold", color="red")

        # chiffres quantiques
        for k in range(1, cq_max + 1):
            plt.text(0.2, h_quantique + 0.3, "n = " + str(k),
                     fontsize=9, weight="bold", style="italic", color="black")
            h_quantique -= 2

        # construction des cases
        for n in range(len(cc)-2):
            for i, j in enumerate(cc):

                # case pour n = 1
                h_quantique = 2 * cq_max - 1
                if cc[i-1] == 1:
                    ax.add_patch(Rectangle((4, h_quantique), 1,
                                 1, fill=False, color="green"))

                    if cc[i-1] == 1 and cc[i] in mSc.keys():
                        if cc[i+1] % 2 == 0:
                            plt.arrow(4.3, h_quantique + 0.1, 0, 0.6,
                                      width=0.035, color="black")
                            plt.arrow(4.65, h_quantique + 0.83, 0, -
                                      0.6, width=0.035, color="black")
                        else:
                            plt.arrow(4.5, h_quantique + 0.1, 0, 0.6,
                                      width=0.035, color="black")
                    plt.text(4.4, h_quantique + 1.2, "s",
                             fontsize=10, color="lightblue")

                # cases pour n = 2
                if cc[i-1] == 2:
                    if cc[i] == "s":
                        ax.add_patch(Rectangle((4, h_quantique - 2),
                                     1, 1, fill=False, color="green"))
                        plt.text(4.4, h_quantique - 0.8, "s",
                                 fontsize=10, color="lightblue")

                        if cc[i+1] % 2 == 0:
                            plt.arrow(4.3, h_quantique - 1.9, 0, 0.6,
                                      width=0.035, color="black")
                            plt.arrow(4.65, h_quantique - 1.17, 0, -
                                      0.6, width=0.035, color="black")
                        else:
                            plt.arrow(4.5, h_quantique - 1.9, 0, 0.6,
                                      width=0.035, color="black")

                    if cc[i] == "p":
                        for j in range(6, 9):
                            ax.add_patch(
                                Rectangle((j, h_quantique - 2), 1, 1, fill=False, color="green"))
                        plt.text(7.4, h_quantique - 0.8, "p",
                                 fontsize=10, color="lightblue")

                        if cc[i+1] // 3 > 0:
                            index = cc[i+1] - 3

                            for i in range(3):
                                plt.arrow(6.3 + i, h_quantique - 1.9,
                                          0, 0.6, width=0.035, color="black")
                            for i in range(index):
                                plt.arrow(6.65 + i, h_quantique - 1.17,
                                          0, -0.6, width=0.035, color="black")

                        else:
                            index = cc[i+1]
                            for i in range(index):
                                plt.arrow(6.3 + i, h_quantique - 1.9,
                                          0, 0.6, width=0.035, color="black")

                # cases pour n = 3
                if cc[i-1] == 3:
                    if cc[i] == "s":
                        ax.add_patch(Rectangle((4, h_quantique - 4),
                                     1, 1, fill=False, color="green"))
                        plt.text(4.3, h_quantique - 2.8, "s",
                                 fontsize=10, color="lightblue")

                        if cc[i+1] % 2 == 0:
                            plt.arrow(4.3, h_quantique - 3.9, 0, 0.6,
                                      width=0.035, color="black")
                            plt.arrow(4.65, h_quantique - 3.17, 0, -
                                      0.6, width=0.035, color="black")
                        else:
                            plt.arrow(4.5, h_quantique - 3.9, 0, 0.6,
                                      width=0.035, color="black")

                    if cc[i] == "p":

                        for j in range(6, 9):
                            ax.add_patch(
                                Rectangle((j, h_quantique - 4), 1, 1, fill=False, color="green"))
                            plt.text(7.3, h_quantique - 2.8, "p",
                                     fontsize=10, color="lightblue")

                        if cc[i+1] // 3 > 0:
                            index = cc[i+1] - 3
                            for i in range(3):
                                plt.arrow(6.3 + i, h_quantique - 3.9,
                                          0, 0.6, width=0.035, color="black")
                            for i in range(index):
                                plt.arrow(6.65 + i, h_quantique - 3.17,
                                          0, -0.6, width=0.035, color="black")
                        else:
                            index = cc[i+1]
                            for i in range(index):
                                plt.arrow(6.3 + i, h_quantique - 3.9,
                                          0, 0.6, width=0.035, color="black")

                    if cc[i] == "d":
                        for j in range(10, 15):
                            ax.add_patch(
                                Rectangle((j, h_quantique - 4), 1, 1, fill=False, color="green"))
                        plt.text(12.3, h_quantique - 2.8, "d",
                                 fontsize=10, color="lightblue")

                        if cc[i+1] // 5 > 0:
                            index = cc[i+1] - 5
                            for i in range(5):
                                plt.arrow(10.3 + i, h_quantique - 3.9,
                                          0, 0.6, width=0.035, color="black")
                            for i in range(index):
                                plt.arrow(10.65 + i, h_quantique - 3.17,
                                          0, -0.6, width=0.035, color="black")
                        else:
                            index = cc[i+1]
                            for i in range(index):
                                plt.arrow(10.3 + i, h_quantique - 3.9,
                                          0, 0.6, width=0.035, color="black")

                # cases pour n = 4
                if cc[i-1] == 4:
                    if cc[i] == "s":
                        ax.add_patch(Rectangle((4, h_quantique - 6),
                                     1, 1, fill=False, color="green"))
                        plt.text(4.3, h_quantique - 4.8, "s",
                                 fontsize=10, color="lightblue")

                        if cc[i+1] % 2 == 0:
                            plt.arrow(4.3, h_quantique - 5.9, 0, 0.6,
                                      width=0.035, color="black")
                            plt.arrow(4.65, h_quantique - 5.17, 0, -
                                      0.6, width=0.035, color="black")
                        else:
                            plt.arrow(4.5, h_quantique - 5.9, 0, 0.6,
                                      width=0.035, color="black")

                    if cc[i] == "p":
                        for j in range(6, 9):
                            ax.add_patch(
                                Rectangle((j, h_quantique - 6), 1, 1, fill=False, color="green"))
                            plt.text(7.3, h_quantique - 4.8, "p",
                                     fontsize=10, color="lightblue")

                        if cc[i+1] // 3 > 0:
                            index = cc[i+1] - 3
                            for i in range(3):
                                plt.arrow(6.3 + i, h_quantique - 5.9,
                                          0, 0.6, width=0.035, color="black")
                            for i in range(index):
                                plt.arrow(6.65 + i, h_quantique - 5.17,
                                          0, -0.6, width=0.035, color="black")
                        else:
                            index = cc[i+1]
                            for i in range(index):
                                plt.arrow(6.3 + i, h_quantique - 5.9,
                                          0, 0.6, width=0.035, color="black")

                    if cc[i] == "d":
                        for j in range(10, 15):
                            ax.add_patch(
                                Rectangle((j, h_quantique - 6), 1, 1, fill=False, color="green"))
                        plt.text(12.3, h_quantique - 2.8, "d",
                                 fontsize=10, color="lightblue")

                        if cc[i+1] // 5 > 0:
                            index = cc[i+1] - 5
                            for i in range(5):
                                plt.arrow(10.3 + i, h_quantique - 5.9,
                                          0, 0.6, width=0.035, color="black")
                            for i in range(index):
                                plt.arrow(10.65 + i, h_quantique - 5.17,
                                          0, -0.6, width=0.035, color="black")
                        else:
                            index = cc[i+1]
                            for i in range(index):
                                plt.arrow(10.3 + i, h_quantique - 5.9,
                                          0, 0.6, width=0.035, color="black")

                    if cc[i] == "f":
                        for j in range(16, 23):
                            ax.add_patch(
                                Rectangle((j, h_quantique - 6), 1, 1, fill=False, color="green"))
                        plt.text(19.4, h_quantique - 4.8, "f",
                                 fontsize=10, color="lightblue")

                        if cc[i+1] // 7 > 0:
                            index = cc[i+1] - 7
                            for i in range(7):
                                plt.arrow(16.3 + i, h_quantique - 5.9,
                                          0, 0.6, width=0.035, color="black")
                            for i in range(index):
                                plt.arrow(16.65 + i, h_quantique - 5.17,
                                          0, -0.6, width=0.035, color="black")
                        else:
                            nombre_ite = cc[i+1]
                            for i in range(index):
                                plt.arrow(16.3 + i, h_quantique - 5.9,
                                          0, 0.6, width=0.035, color="black")

                # cases pour n = 5
                if cc[i-1] == 5:
                    if cc[i] == "s":
                        ax.add_patch(Rectangle((4, h_quantique - 8),
                                     1, 1, fill=False, color="green"))
                        plt.text(4.3, h_quantique - 6.8, "s",
                                 fontsize=10, color="lightblue")

                        if cc[i+1] % 2 == 0:
                            plt.arrow(4.3, h_quantique - 7.9, 0, 0.6,
                                      width=0.035, color="black")
                            plt.arrow(4.65, h_quantique - 7.17, 0, -
                                      0.6, width=0.035, color="black")
                        else:
                            plt.arrow(4.5, h_quantique - 7.9, 0, 0.6,
                                      width=0.035, color="black")

                    if cc[i] == "p":
                        for j in range(6, 9):
                            ax.add_patch(
                                Rectangle((j, h_quantique - 8), 1, 1, fill=False, color="green"))
                            plt.text(7.3, h_quantique - 6.8, "p",
                                     fontsize=10, color="lightblue")

                        if cc[i+1] // 3 > 0:
                            index = cc[i+1] - 3
                            for i in range(3):
                                plt.arrow(6.3 + i, h_quantique - 7.9,
                                          0, 0.6, width=0.035, color="black")
                            for i in range(index):
                                plt.arrow(6.65 + i, h_quantique - 7.17,
                                          0, -0.6, width=0.035, color="black")
                        else:
                            index = cc[i+1]
                            for i in range(index):
                                plt.arrow(6.3 + i, h_quantique - 7.9,
                                          0, 0.6, width=0.035, color="black")

                    if cc[i] == "d":
                        for j in range(10, 15):
                            ax.add_patch(
                                Rectangle((j, h_quantique - 8), 1, 1, fill=False, color="green"))
                        plt.text(12.3, h_quantique - 6.8, "d",
                                 fontsize=10, color="lightblue")

                        if cc[i+1] // 5 > 0:
                            index = cc[i+1] - 5
                            for i in range(5):
                                plt.arrow(10.3 + i, h_quantique - 7.9,
                                          0, 0.6, width=0.035, color="black")
                            for i in range(index):
                                plt.arrow(10.65 + i, h_quantique - 7.17,
                                          0, -0.6, width=0.035, color="black")
                        else:
                            index = cc[i+1]
                            for i in range(index):
                                plt.arrow(10.3 + i, h_quantique - 7.9,
                                          0, 0.6, width=0.035, color="black")

                    if cc[i] == "f":
                        for j in range(16, 23):
                            ax.add_patch(
                                Rectangle((j, h_quantique - 8), 1, 1, fill=False, color="green"))
                        plt.text(19.4, h_quantique - 4.8, "f",
                                 fontsize=10, color="lightblue")

                        if cc[i+1] // 7 > 0:
                            index = cc[i+1] - 7
                            for i in range(7):
                                plt.arrow(16.3 + i, h_quantique - 7.9,
                                          0, 0.6, width=0.035, color="black")
                            for i in range(index):
                                plt.arrow(16.65 + i, h_quantique - 7.17,
                                          0, -0.6, width=0.035, color="black")
                        else:
                            nombre_ite = cc[i+1]
                            for i in range(index):
                                plt.arrow(16.3 + i, h_quantique - 7.9,
                                          0, 0.6, width=0.035, color="black")

                # cases pour n = 6
                if cc[i-1] == 6:
                    if cc[i] == "s":
                        ax.add_patch(Rectangle((4, h_quantique - 10),
                                     1, 1, fill=False, color="green"))
                        plt.text(4.3, h_quantique - 8.8, "s",
                                 fontsize=10, color="lightblue")

                        if cc[i+1] % 2 == 0:
                            plt.arrow(4.3, h_quantique - 9.9, 0, 0.6,
                                      width=0.035, color="black")
                            plt.arrow(4.65, h_quantique - 9.17, 0, -
                                      0.6, width=0.035, color="black")
                        else:
                            plt.arrow(4.5, h_quantique - 9.9, 0, 0.6,
                                      width=0.035, color="black")

                    if cc[i] == "p":
                        for j in range(6, 9):
                            ax.add_patch(
                                Rectangle((j, h_quantique - 10), 1, 1, fill=False, color="green"))
                            plt.text(7.3, h_quantique - 8.8, "p",
                                     fontsize=10, color="lightblue")

                        if cc[i+1] // 3 > 0:
                            index = cc[i+1] - 3
                            for i in range(3):
                                plt.arrow(6.3 + i, h_quantique - 9.9,
                                          0, 0.6, width=0.035, color="black")
                            for i in range(index):
                                plt.arrow(6.65 + i, h_quantique - 9.17,
                                          0, -0.6, width=0.035, color="black")
                        else:
                            index = cc[i+1]
                            for i in range(index):
                                plt.arrow(6.3 + i, h_quantique - 9.9,
                                          0, 0.6, width=0.035, color="black")

                    if cc[i] == "d":
                        for j in range(10, 15):
                            ax.add_patch(
                                Rectangle((j, h_quantique - 10), 1, 1, fill=False, color="green"))
                        plt.text(12.3, h_quantique - 8.8, "d",
                                 fontsize=10, color="lightblue")

                        if cc[i+1] // 5 > 0:
                            index = cc[i+1] - 5
                            for i in range(5):
                                plt.arrow(10.3 + i, h_quantique - 9.9,
                                          0, 0.6, width=0.035, color="black")
                            for i in range(index):
                                plt.arrow(10.65 + i, h_quantique - 9.17,
                                          0, -0.6, width=0.035, color="black")
                        else:
                            index = cc[i+1]
                        for i in range(index):
                            plt.arrow(10.3 + i, h_quantique - 9.9, 0,
                                      0.6, width=0.035, color="black")

                # cases pour n = 7
                if cc[i-1] == 7:
                    if cc[i] == "s":
                        ax.add_patch(Rectangle((4, h_quantique - 12),
                                     1, 1, fill=False, color="green"))
                        plt.text(4.4, h_quantique - 10.8, "s",
                                 fontsize=10, color="lightblue")

                        if cc[i+1] % 2 == 0:
                            plt.arrow(4.3, h_quantique - 11.9, 0,
                                      0.6, width=0.035, color="black")
                            plt.arrow(4.65, h_quantique - 11.17,
                                      0, -0.6, width=0.035, color="black")
                        else:
                            plt.arrow(4.5, h_quantique - 11.9, 0,
                                      0.6, width=0.035, color="black")

                    if cc[i] == "p":
                        for j in range(6, 9):
                            ax.add_patch(
                                Rectangle((j, h_quantique - 12), 1, 1, fill=False, color="green"))
                        plt.text(7.4, h_quantique - 10.8, "p",
                                 fontsize=10, color="lightblue")

                        if cc[i+1] // 3 > 0:
                            index = cc[i+1] - 3

                            for i in range(3):
                                plt.arrow(6.3 + i, h_quantique - 11.9,
                                          0, 0.6, width=0.035, color="black")
                            for i in range(index):
                                plt.arrow(6.65 + i, h_quantique - 11.17,
                                          0, -0.6, width=0.035, color="black")

                        else:
                            index = cc[i+1]
                            for i in range(index):
                                plt.arrow(6.3 + i, h_quantique - 11.9,
                                          0, 0.6, width=0.035, color="black")

        plt.axis('equal')
        plt.axis('off')
        plt.show()

    # Appel des fonctions

    if action == "structureComplete":
        return structureComplete(X)[0]

    elif action == "structureAllegee":
        return structureAllegee(X)

    elif action == "configurationValence":
        return configurationValence(X)[0]

    elif action == "casesQuantiques":
        return casesQuantiques(X)

# Appel de la fonction atomistique


print(atomistique("Mg", "structureComplete"))
