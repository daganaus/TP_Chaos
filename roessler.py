# -*- coding: utf-8 -*-
# Ce code est destiné au TP numérique de Chaos en M1 Physique
# de l'Université Grenoble Alpes. 
# Il a été écrit par Vincent Rossetto (vincent.rossetto@grenoble.cnrs.fr)

## INSTRUCTIONS ##
# Pour le travail de TP demandé, ce code n'est qu'un point de départ.
# Il permet de réaliser des opérations simples mais tel quel, il ne permet
# pas de répondre aux questions. Pour répondre aux questions, il faut le
# modifier en changeant des paramètres ou en ajoutant des fonctions.
# 
# Le programme obtenu pour chaque réponse devra se baser sur celui-ci,
# ou sur les programmes des questions précédentes. À chaque 
# question, faire une sauvegarde du code utilisé, en le nommant par exemple
#   TP_chaos-1.py 
# pour le programme qui répond à la question 1.
# Faites des sauvegardes régulièrement pour éviter de perdre le travail
# effectué. 

# Des documents sur python sont disponibles aux adresses suivantes :
# Sur le site de l'UE Chaos et applications :
#    
# Un petit cours rapide d'introduction / rappels sur le langage :    
#   https://cours.univ-grenoble-alpes.fr/mod/resource/view.php?id=505241
#
# Un cours de calcul scientifique avec python, niveau L3 physique
#   https://cours.univ-grenoble-alpes.fr/mod/resource/view.php?id=505243
#
# Le site de référence en français
#   https://docs.python.org/fr/3
#
# Le code ci-dessous contient beaucoup de commentaires destinés à vous aider
# à comprendre comment il fonctionne, quelles sont les fonctions # utilisées. 
# N'écrivez pas de commentaires aussi détaillés pour votre code ; cependant 
#
############################################################################## 
#           IL EST INDISPENSABLE DE COMMENTER VOTRE CODE !                   #
##############################################################################
#
# Bon TP !
# Vincent Rossetto
# 2025-04-11

# Importation des paquets standards pour le calcul
# Numpy est renommé np, comme usuellement
import  numpy as np

# Importation de la fonction odeint de scipy qui permet d'intégrer 
# l'équation numériquement
from scipy.integrate import odeint

# Les outils d'affichage de matplotlib pour python
# Le paquet est renommé plt
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d 
# Mise en activité du mode interactif
plt.ion()
# Définition de la méthode d'affichage des graphiques.
# Entrer '%matplotlib qt' dans la console si la fenêtre ne s'ouvre pas.
plt.switch_backend('Qt5Agg')

# Importation de quelques widgets pour manipuler les graphiques et
# les images. Voir plus bas comment les utiliser (c'est facile !)
from matplotlib.widgets import Slider, Button

# DÉFINITION DE LA FIGURE
# On commence par ouvrir une figure à l'aide du module plt
# qui est matplotlib.pyplot
fig=plt.figure()
# Cette fonction récupère la définition des axes de la figure.
ax=fig.add_subplot(projection='3d')
# Avec une ancienne version de matplotlib, remplacer la ligne ci-dessus par
# ax=fig.gca(projection='3d')
 
# PARAMÈTRES 
# Nombre de points calculés
N = 10000
# Définition des paramètres du système dynamique
# Paramètres initiaux (a,b,c) 
(a,b,c)=(0.25, 1 , 2.4) 
# a et b gardent ces valeurs, mais c sera modifié.
c0=c
# Temps initial du tracé
t0=0
# Position initiale (x,y,z)
R_in = [1, 0.5, 0.5]   # ATTENTION : condition initiale "mauvaise"
# Durée de la solution tracée
t1=100

# La fonction qui définit le système dynamique:
# Il s'agit simplement qui à un vecteur R=(x,y,z) associe
# les dérivées de (x,y,z). 
def Roessler(R, t, a, b, c):
  """Les équations de Rössler"""
  return [-R[1] - R[2],
          R[0] + a*R[1],
          b + (R[0] -c)*R[2]]

# Une fonction auxiliaire qui donne la solution
# des équations à partir d'une condition initiale
# pour une certaine durée. 
def solve_Roessler(r0, duree, parametres) :
  """La solution des équations de Rössler, calculée
     pour une durée donnée à partir de la condition
     initiale r0.
     Le résultat contient un vecteur temps et un
     tableau des positions."""
  # Définition des temps
  t=np.linspace(0, duree, N)
  # Calcul de la solution jusqu'au temps demandé .
  # Le résultat est un tableau de valeurs R[0], ... , R[N]
  # tels que R[i] est le vecteur position au temps i.
  # R[i,j] est la coordonnée j du vecteur R[i].
  R=odeint(Roessler, r0, t, args=parametres)
  # La fonction retourne les temps et les positions.
  return t, R

def Poincare(R, parametres) :
  """Fonction de tracé de la section de Poincaré
  pour le système de Rössler"""
  # On ne garde que les points où la coordonnée x est nulle.
  # On ne garde que les points où la coordonnée y est positive
  # et on ne garde que les points où la coordonnée z est positive
  # (ce qui correspond à une section de Poincaré).
  # On utilise le fait que R est un tableau numpy.
  # On peut donc faire des opérations sur les tableaux.
  R = R[(R[:,0] < 0.01) & (R[:,0] > -0.01) & (R[:,1] > 0) & (R[:,2] > 0)]
  # On trace la section de Poincaré
  ax.plot3D(R[:,1], R[:,2], R[:,0], 'green')
  # On ajoute le point fixe sur la figure.
  Roessler_fixed_point(parametres)

# Le point fixe du système de Rössler
# Il correspond au point tel que Roessler(R,t,a,b,c)=0
def Roessler_fixed_point(parametres) :
  a=parametres[0]
  c=parametres[2]
  D = c**2-4*a*parametres[1]
  xp0 = (c - np.sqrt(D))/2
  xp1 = -xp0/a
  xp2 = xp0/a
  
  # Tracé du point fixe en rouge
  ax.plot3D([xp0], [xp1], [xp2],  marker='.', linestyle='none', color='red')

def trace_Roessler(r0,t0,t1,parametres) :
  """Fonction de tracé d'une trajectoire du système de Rössler
  avec les paramètres a, b et c, pour une durée t1 à partir du temps t0
  et de la position initiale R0"""
  # On avance dans le temps d'une durée t0.
  # Le signe _ signifie qu'on ne garde pas en mémoire
  # le vecteur temps (il est inutile pour tracer cette courbe.
  _, r=solve_Roessler(r0,t0,parametres)
  # On récupère la dernière valeur calculée qui sera la première
  # valeur tracée.
  r1=r[-1]
  # On résoud de nouveau à partir de la nouvelle origine.
  _, R=solve_Roessler(r1,t1,parametres)
  # R.T est la transposée de R, donc R.T[j,i] est la coordonnée j au temps i.
  # Ainsi R.T[0] est le vecteur des valeur de la coordonnée 0 (soit x)
  # à tous les temps i.
  [X, Y, Z] = R.T 
  # On efface la figure avant de tracer.
  ax.clear()
  # On trace la courbe calculée.
  ax.plot3D(X, Y, Z, 'blue')
  # On ajoute le point fixe sur la figure.
  Roessler_fixed_point(parametres)

# FONCTIONS POUR LES WIDGETS
# La fonction quitter ne fait que fermer la fenêtre en cours d'utilisation.
# Cette fonction est nécessaire pour créer un bouton qui effectue 
# l'action de fermer la fenêtre.
def quitter(_):
  plt.close()

# Fonction de mise à jour de l'affichage, pour prendre compte la modification 
# d'un paramètre. Elle sera activée à chaque modification du paramètre c
# ou du temps t0 (voir plus bas).
def update(_):
  # On récupère la valeur du paramètre c indiqué par la barre de glissement.
  c=barre_c.val
  # On récupère la valeur du temps indiqué par la barre de glissement.
  t0=barre_t.val
  # On trace la nouvelle figure
  trace_Roessler(R_in,t0,t1,(a,b,c))

# Fonction de remise à zéro des glisseurs.
# Une fois cela fait, on réactualise l'affichage
def reset(_):
  barre_c.reset()
  barre_t.reset()
  update(0)

# TRACÉ DES WIDGETS
# Dessin de la barre de glissement pour le paramètre c
# Les nombres entrés sont les coordonnées du rectangle contenant la barre.
axe_c = plt.axes([0.1, 0.06, 0.65, 0.03])
# On crée ensuite un widget de type slider avec le rectangle défini axe_c.
# Les valeurs de la barre vont de 1 (à gauche) à 15 (à droite).
# Le nom indiqué à gauche de la barre est 'c' et sa valeur initiale est c0.
barre_c= Slider(axe_c, 'c', 1 , 5, c0)

# Dessin de la barre de glissement pour le paramètre t0
axe_t = plt.axes([0.1, 0.02, 0.65, 0.03])
# Widget de type slider pour le paramètre t0 (la valeur initiale est t0=0).
barre_t= Slider(axe_t,'t0',0,10000)

# Dessin d'un rectangle qui sera la bouton action.
cadre_raz=plt.axes([0.85, 0.06, 0.1, 0.03])
# Widget de type bouton 
bouton_raz=Button(cadre_raz,'R. à 0')

# Dessin d'un rectangle qui sera la bouton Fin
cadre_fin = plt.axes([0.85, 0.02, 0.1, 0.03])
# Widget de type bouton 
bouton_fin=Button(cadre_fin,'Fin')

# Les widgets sont créés, mais il faut maintenant associer 
# ce qu'il se passe quand on les utilise.

# ACTIVATION DES WIDGETS 
# Lorsque l'on change la valeur de la barre (on_changed) l'action update
# est exécutée.
barre_c.on_changed(update)
# Même chose pour la barre de t0
barre_t.on_changed(update)
# Lorsque l'on clique sur remise à zéro, on réinitialise
bouton_raz.on_clicked(reset)
# Lorsque l'on clique sur Fin, on ferme la fenêtre
bouton_fin.on_clicked(quitter)

# Initialisation du programme
# On fait une mise à jour avec update et un paramètre (n'importe lequel).
# Cela permet de mettre les valeurs des paramètres a,b,c,t0 à jour,
# en accord avec les valeurs indiquées sur les glisseurs.
update(0)

# AFFICHAGE DE LA FENÊTRE
plt.show(block=True)
# Le programme "tourne" maintenant dans la fenêtre. Il attend
# que les widgets soient utilisés pour effectuer les actions
# demandées. Tant que la fenêtre est ouverte, le programme continue 
# d'attendre et exécuter les opérations controlées par les widgets.
# Si on ferme la fenêtre, comme il n'y a plus rien en dessous de ce code,
# python a réalisé toutes les opérations demandées et termine le programme.
