import  numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d 
plt.ion()
plt.switch_backend('Qt5Agg')
from matplotlib.widgets import Slider, Button

fig=plt.figure(figsize=(8,6))
ax=fig.add_subplot(projection='3d')
# Avec une ancienne version de matplotlib, remplacer la ligne ci-dessus par
#ax=fig.gca(projection='3d')
ax.set_position([0, 0.15, 1, 0.8])

fig2, ax2D = plt.subplots()
 
# PARAMÈTRES 
# Nombre de points calculés
N = 10000
# Définition des paramètres du système dynamique
# Paramètres initiaux (a,b,c) 
(a,b,c)=(0.25, 1, 2.4) 
# a et b gardent ces valeurs, mais c sera modifié.
c0=c
# Temps initial du tracé
t0=0
# Position initiale (x,y,z)
R_in = [0, 1, 0.5]   # ATTENTION : condition initiale "mauvaise"
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

## Calcul de la section de Point carré

# équation du plan (section de Poincaré) : ax + by + cz + d = 0
# n : vecteur normal au plan, n = (a,b,c)
# d : position du plan
# r1, r2 : deux points successifs de la trajectoire
def test_intersect(r1, r2, n, d): # fonction qui teste si le segment [r1,r2] intersecte le plan
  if np.dot(r1,n)*np.dot(r2,n) < 0 : # On regarde si r1.n et r2.n ont des signes opposés
      return True
  else:
      return False

def intersect(r1, r2, n, d): # fonction qui retourne le point d'intersection s'il existe
    direction = r2-r1
    denom = np.dot(direction, n)
    if test_intersect(r1, r2, n, d): #Il existe une intersection on la calcule en résolvant n.dir = d
        t = -(np.dot(r1,n) + d)/denom
        if 0 <= t <= 1:
            return r1 + t*direction
    
def Poincare(R, n, d):
    P_intersect = [] # liste des points d'intersection de la trajectoire avec la section de Poincaré
    for i in range(len(R)-1):
        pt = intersect(R[i], R[i+1], n, d)
        if pt is not None:
            P_intersect.append(pt)
    #print(np.array(P_intersect))
    return np.array(P_intersect)



        
    

# Une fonction auxiliaire qui donne la solution
# des équations à partir d'une condition initiale
# pour une certaine durée. 
def solve_Roessler(r0, parametres, duree, npoints=N) :
  """La solution des équations de Rössler, calculée
     pour une durée donnée à partir de la condition
     initiale r0.
     Le résultat contient un vecteur temps et un
     tableau des positions."""
  # Définition des temps
  t=np.linspace(0, duree, npoints)
  R=odeint(Roessler, r0, t, args=parametres)
  # La fonction retourne les temps et les positions.
  return t, R

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

def trace_Roessler(r0, parametres, t0, t1, npoints=N) :
  """Fonction de tracé d'une trajectoire du système de Rössler
  avec les paramètres a, b et c, pour une durée t1 à partir du temps t0
  et de la position initiale R0"""
  # On avance dans le temps d'une durée t0
  # le nombre de pas est réduit proportionnellement à t0 pour éviter
  # de prendre trop de temps sur cette initialisation.
  n=int(t0/t1*npoints) + 1
  # Le signe _ signifie qu'on ne garde pas en mémoire
  # le vecteur temps (il est inutile pour tracer cette courbe)
  _, r=solve_Roessler(r0, parametres, t0, n)
  # On récupère la dernière valeur calculée qui sera la première
  # valeur tracée.
  r1=r[-1]
  # On résoud de nouveau à partir de la nouvelle origine.
  _, R=solve_Roessler(r1, parametres, t1, npoints)
  P_intersect = Poincare(R, [0,1,0], -1)
  X_intersect, Y_intersect, Z_intersect = P_intersect.T
  ax2D.clear()
  ax2D.plot(X_intersect, Z_intersect, '.', color='darkorange', label="Section de Poincaré (y=0)")
  ax2D.set_xlabel("x")
  ax2D.set_ylabel("z")
  ax2D.set_title("Section de Poincaré 2D")
  ax2D.legend()
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
  t0=barre_t0.val
  t1=barre_t1.val
  # On récupère le nombre de pas 
  n=np.power(10., barre_N.val)
  barre_N.valtext.set_text(f"{n:5.1e}")
  # On trace la nouvelle figure
  trace_Roessler(R_in, (a,b,c), t0, t1, int(n))

# Fonction de remise à zéro des glisseurs.
# Une fois cela fait, on réactualise l'affichage
def reset(_):
  barre_c.reset()
  barre_t0.reset()
  barre_t1.reset()
  barre_N.reset()
  update(0)

# TRACÉ DES WIDGETS
# Dessin de la barre de glissement pour le paramètre c 
# Les nombres entrés sont les coordonnées du rectangle contenant la barre.
axe_c = plt.axes([0.1, 0.07, 0.65, 0.04])
# On crée ensuite un widget de type slider avec le rectangle défini axe_c.
# Les valeurs de la barre vont de 1 (à gauche) à 15 (à droite).
# Le nom indiqué à gauche de la barre est 'c' et sa valeur initiale est c0.
barre_c= Slider(ax=axe_c, label='c', valmin=2, valmax=5, valinit=c0, track_color='darkgreen')

# Dessin de la barre de glissement pour les paramètres de temps (t0 et t1)
axe_t0 = plt.axes([0.1, 0.01, 0.25, 0.03])
axe_t1 = plt.axes([0.45, 0.01, 0.30, 0.03])
# Widget de type slider pour les temps (la valeur initiale est t0=0).
barre_t0= Slider(ax=axe_t0, label='t0', valmin=0, valmax=10000, valinit=0)
barre_t1= Slider(ax=axe_t1, label='t1', valmin=0, valmax=10000, valinit=100)

# De même avec le paramètre N (il est en échelle logarithmique de base 10)
axe_N = plt.axes([0.1, 0.04, 0.65, 0.03])
barre_N= Slider(ax=axe_N, label='N', valmin=2, valmax=6, valinit=4)  
# Attention à ne pas choisir N trop grand !

# Dessin d'un rectangle qui sera la bouton action.
cadre_raz=plt.axes([0.85, 0.05, 0.1, 0.03])
# Widget de type bouton 
bouton_raz=Button(cadre_raz,'R. à 0')

# Dessin d'un rectangle qui sera la bouton Fin
cadre_fin = plt.axes([0.85, 0.01, 0.1, 0.03])
# Widget de type bouton 
bouton_fin=Button(cadre_fin,'Fin')

barre_c.on_changed(update)
barre_t0.on_changed(update)
barre_t1.on_changed(update)
barre_N.on_changed(update)
bouton_raz.on_clicked(reset)
bouton_fin.on_clicked(quitter)
update(0)

plt.show(block=True)



