import numpy as np
from math import *
import os
import linecache
import sys
GREAT=1e10
eps=1e-5
#___________________________________________________________Parametres utilisateur__________________________________________________#
# Version de castem (a modifier)
Version_castem=16
# Pas de temps pour le calcul thermo-mecanique
Pas=1 # seconde                                                                                                                                                                    A____ B
run=False # Lancement automatique du calcul dans Cast3m a la fin de son ecriture. False est conseille                                                                               D1\ | 
MaillageCartesien=True # Si le maillage est cartesien, des coins peuvent apparaitre et doivent etre supprimes (divergence des calculs mecaniques). Dans ce cas, on casse les angles    \|D2
rupelian=True # Choix entre deux materiaux lors du calcul thermo-mecanique.                                                    #                                                        | 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                    Parametres du maillage utilisateurs                       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Nombre de mailles par ligne
Discretisation_Surf_OpenFOAM=150 # Discretisation fine de la zone d'interet
Discretisation_sides=20 # Discretisation grossiere d'une zone sans interet
Discretisation_haut=50 # Discretisation grossiere d'une zone sans interet
# Hauteur d'extrusion
H=2.4
# Direction d'extrusion
Dir_ex=np.zeros(3)
Dir_ex[0],Dir_ex[1],Dir_ex[2]=0,0,1
# Vecteur normal au plan de coupe
vect_normal=np.zeros(3) 
vect_normal[0],vect_normal[1],vect_normal[2]=0,1,0
#Point A appartenant au plan de coupe
Pt_plan=np.zeros(3)
Pt_plan[0],Pt_plan[1],Pt_plan[2]=1252,1285.95,192
# Altitude dans la direction Dir_ex en-dessous de laquelle la geometrie n'est pas gardee
CoupeMini=185.8
# Altitude dans la direction Dir_ex au-delà de laquelle la geometrie n'est pas gardee
CoupeMax=190

#___________________________________________________________________________________________________________________________________#

# Produit vectoriel entre Dir_ex et vect_normal: c'est donc le deuxieme vecteur du plan avec Dir_ex. (Dir_ortho,Dir_ex) forme une base orthonormale du plan de coupe
Dir_ortho=np.cross(Dir_ex,vect_normal)
# Matrice de changement de base
P=np.zeros((3,3))
for i in range(0,3):
  P[0,i]=Dir_ortho[i]
  P[1,i]=Dir_ex[i]
  P[2,i]=vect_normal[i]
  
# Fonction produit matriciel entre une matrice et un vecteur
def ProdMat(A,B):
  eps=1e-3
  C=np.zeros(3)
  for i in range(0,3):
    for j in range(0,3):
      C[i]+=A[i,j]*B[j]
    if abs(C[i])<eps:
      C[i]=0
  return C

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                        Definition du plan de coupe                           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Equation cartesienne du plan passant par A, le point appartenant au plan de coupe
def Eq_Plan(x,y,z):
  a=vect_normal[0]
  b=vect_normal[1]
  c=vect_normal[2]
  d=-(a*Pt_plan[0]+b*Pt_plan[1]+c*Pt_plan[2])
  return a*x+b*y+c*z+d
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#           Fonction renvoyant vrai si l'element est dans la liste             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def Is_In(list,point):
  len_list=len(list)
  for i in range(0,len_list):
    if (point[0]-list[i][0])**2+(point[1]-list[i][1])**2+(point[2]-list[i][2])**2<eps:
      return True
  return False

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  Obtention du numero de la premiere face et du nombre de faces des surfaces  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
f = open('constant/polyMesh/boundary','r')
contenu = f.read()

Name_surf=[] # Stockage des noms de toutes les surfaces
Shouldbreak=False
for k in range(0,len(contenu)):
  if contenu[k]=="(":
    break
i=1
while 1!=2:
  while (contenu[k+i]!="{"):
    i+=1
    if k+i>=len(contenu)-1: # Si on arrive a la fin du fichier, on arrete
      Shouldbreak=True
      break
  if Shouldbreak==True:
    break
  Name_surf1=[]
  l=0
  while contenu[k+i-6-l]!=" ":
    Name_surf1+=contenu[k+i-6-l]
    l+=1
  Name_surf1=str("".join(reversed(Name_surf1)))
  Name_surf.append(Name_surf1) # Contient tous les noms de surface
  i+=1

# Recuperation du nombre de faces associees a chaque surface ainsi que le numero des 1eres faces
nb_surf=len(Name_surf)

nFaces=np.zeros(nb_surf)
startFace=np.zeros(nb_surf)
for i in range(0,nb_surf):
  Nb_surf=contenu.find(Name_surf[i])
  contenu_2=contenu[Nb_surf:]
  Nb_faces=contenu_2.find('nFaces')
  contenu_3=contenu_2[Nb_faces:]
  Nb_point=contenu_3.find(';')
  contenu_4=contenu_3[Nb_point+1:]
  Nb_point_2=contenu_4.find(';')
  ncar=0
  nbFaces=""
  car=""
  while car!=" ":
    ncar+=1
    car=contenu_3[Nb_point-ncar]
    nbFaces+=car
  car=""
  ncar=0
  nstartFace=""
  while car!=" ":
    ncar+=1
    car=contenu_4[Nb_point_2-ncar]
    nstartFace+=car    
  nFaces[i]=int("".join(reversed(nbFaces))) # Nombre de faces de chaque surface i
  startFace[i]=int("".join(reversed(nstartFace))) # Numero de la 1ere face de toutes les surfaces i


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Reecriture du fichier faces dans un fichier faces2 sans superflu       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def Create():
	k=0
	length=0
	Fich='faces2'
	fichier2=open(Fich,'w')
	fichier=open('constant/polyMesh/faces','r')
	lignes=fichier.readlines()
	for ligne in lignes:
		length=length+1
	for ligne in lignes:
		k+=1
		ligne=ligne.replace('(',' ')
		ligne=ligne.replace(')','')
		if k>20 and k<length:
			fichier2.write(ligne)
	fichier2.close
	fichier.close()

Create()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Determine la 1ere face a selectionner dans le fichier faces2           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
Table=np.fromfile("faces2",sep=" ")

ligne=0
k=0
i=0
FirstFace=np.zeros(nb_surf)
for j in range(0,nb_surf):
  while 1!=2:
	  i=i+k
	  if Table[i]>=3 and Table[i]<=8:
		  k=Table[i]+1
	  ligne+=1
	  if ligne==startFace[j]+1:
		  FirstFace[j]=i # Premiere face
		  break


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#               Attribution de chaque numero de face a ListFace                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
ListFace=[]
for j in range(0,nb_surf):
  ListFacePoub=np.zeros((nFaces[j],9))
  l=FirstFace[j]+1
  ListFacePoub[0,0]=Table[l]
  ligne=0
  k=0
  i=0
  er=Table[l-1]
  while ligne<nFaces[j]:
    if Table[l]>=3 and Table[l]<=9 and er==0:
		  ligne+=1
		  k=0
		  i+=1
		  er=Table[l]
    else:
		  ListFacePoub[i,k]=Table[l]
		  k+=1
		  er-=1
    l+=1
    if l==len(Table)-1:
      break
  ListFace.append(ListFacePoub)
os.remove("faces2")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Reecriture du fichier points dans un fichier points2 sans superflu     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def Create2():
	k=0
	length=0
	Fich='points2'
	fichier2=open(Fich,'w')
	fichier=open('constant/polyMesh/points','r')
	lignes=fichier.readlines()
	for ligne in lignes:
		length=length+1
	for ligne in lignes:
		k=k+1
		ligne=ligne.replace('(','')
		ligne=ligne.replace(')','')
		if k>20 and k<length:
			fichier2.write(ligne)
	fichier2.close
	fichier.close()

Create2()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                   Ecriture du fichier Maillage.txt                           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
Maillage='Maillage.txt' # Contient les informations sur les mailles composant les surfaces d'OpenFOAM
fichier=open(Maillage,'w')
for j in range(0,nb_surf):
  for i in range(0,int(nFaces[j])):
	  b=""
	  for k in range(0,10):
		  if int(ListFace[j][i,k])!=0:
			  a=linecache.getline("points2",int(ListFace[j][i,k])+1)		
			  b="{0} {1}".format(i,a)
			  fichier.write(b)
		  else:
			  break

os.remove("points2")
fichier.close()
Table3=np.fromfile("Maillage.txt",count=-1,sep=" ")
os.remove("Maillage.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                   Ecriture du fichier Maillages.geo                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
Gmsh='CASTEM/Maillages.geo'
fichier=open(Gmsh,'w') # Fichier d'entree de GMSH (en l'absence de lissage)
Gmsh_points='Point.txt'
fichier2=open(Gmsh_points,'w') # Collecte les points d'intersection entre le plan de coupe et la surface
Castem_faces='Castem.txt'
fichier3=open(Castem_faces,'w')
n=np.zeros(8)
NB=0
Nb_precedent=0 # Correspond au numero de la face actuellement testee

# Determination de la nature (triangle, carré, etc) de la derniere maille pour arreter la boucle ci-apres
Last=1
while Table3[len(Table3)-4]==Table3[len(Table3)-4-4*Last]:
  Last+=1
# Boucle for qui determine les points a l'intersection de la surface d'OpenFOAM et du plan de coupe
cont=[]
contenu=[]
number=0
while Nb_precedent<len(Table3)-4*Last:
  j=0
  if Table3[Nb_precedent]==0:
    num_face=0
    number+=1 # Correspond au numero de la surface
  while Table3[Nb_precedent+4*j]==Table3[Nb_precedent+4*j+4]:
    n[j]=Eq_Plan(Table3[Nb_precedent+4*j+1],Table3[Nb_precedent+4*j+2],Table3[Nb_precedent+4*j+3])
    j+=1
    if (i==len(Table3)/16-2 and j==3):
      break
  n[j]=Eq_Plan(Table3[Nb_precedent+4*j+1],Table3[Nb_precedent+4*j+2],Table3[Nb_precedent+4*j+3])
  list_test=[n[k] for k in range(0,j+1)]
  test=min(list_test)*max(list_test)
  if test<0: # Dans ce cas, sur la meme face, il y a au moins un point de chaque cote du plan donc le plan intersecte cette face
    # Stockage du numero de face et du barycentre de la surface dans un fichier texte utilise par le second fichier python
    BaryX=np.mean([Table3[Nb_precedent+4*i+1] for i in range(0,j+1)])
    BaryY=np.mean([Table3[Nb_precedent+4*i+2] for i in range(0,j+1)])                                                                                             # A__\____ B
    BaryZ=np.mean([Table3[Nb_precedent+4*i+3] for i in range(0,j+1)])                                                                                             #  |  \   |
    Bary=ProdMat(P,[BaryX-Pt_plan[0],BaryY-Pt_plan[1],BaryZ-Pt_plan[2]]) # Changement de base --> on se place dans (Dir_ortho,Dir_ex,vect_normal)                 #  |   \  |      On garde les deux points d'intersection                                                                                                                                        #  |____\_|
    if (Bary[1]+np.dot(Pt_plan,Dir_ex)<CoupeMini) or (Bary[1]+np.dot(Pt_plan,Dir_ex)>CoupeMax) or (Bary[0]+np.dot(Pt_plan,Dir_ortho)>-1252.2) or\                 # D        C
    (Bary[0]+np.dot(Pt_plan,Dir_ortho)<-1253.4): # En dehors de la zone a garder                                                                                  
      pass                                                                                                                                                        
    else:
      fichier3.write(str(number)+" "+str(num_face)+" "+str(Bary[0])+" "+str(Bary[1])+" "+str(Bary[2])+"\n")
      #Stockage des longueurs des cotes dans length (permet d eviter de tester l'intersection avec AC et BD dans un carre ABCD)
      length=[]
      for l in range(0,j+1):
        if l<j:
          length.append((Table3[Nb_precedent+4*l+1]-Table3[Nb_precedent+4*l+5])**2+(Table3[Nb_precedent+4*l+2]-Table3[Nb_precedent+4*l+4+2])**2+(Table3[Nb_precedent+4*l+3]-Table3[Nb_precedent+4*l+4+3])**2)
        else:
          length.append((Table3[Nb_precedent+4*l+1]-Table3[Nb_precedent+1])**2+(Table3[Nb_precedent+4*l+2]-Table3[Nb_precedent+2])**2+(Table3[Nb_precedent+4*l+3]-Table3[Nb_precedent+3])**2)
      for k in range(0,j+1):
        if k!=j:                                                                                                                                                                     
          if list_test[k]*list_test[k+1]<0: # Dans ce cas, le plan coupe la maille entre ces deux points                                                                                                                    
            x_a,y_a,z_a=Table3[Nb_precedent+4*k+1], Table3[Nb_precedent+4*k+2], Table3[Nb_precedent+4*k+3]                                                                          
            x_b,y_b,z_b=Table3[Nb_precedent+4*k+4+1], Table3[Nb_precedent+4*k+4+2], Table3[Nb_precedent+4*k+4+3]                          
            norme=(x_a-x_b)**2+(y_a-y_b)**2+(z_a-z_b)**2
            if norme in length: # Dans ce cas il y a bien intersection entre le plan et une arete (ce n'est pas le cas si norme n'est pas dans length car dans ce cas il y a intersection entre le plan et AC par ex pour une cell ABCD
              a,b,c=x_b-x_a, y_b-y_a, z_b-z_a                                                                                                                                                                                 
              t=-(vect_normal[0]*x_a+vect_normal[1]*y_a+vect_normal[2]*z_a-(vect_normal[0]*Pt_plan[0]+vect_normal[1]*Pt_plan[1]+vect_normal[2]*Pt_plan[2]))/(a*vect_normal[0]+b*vect_normal[1]+c*vect_normal[2])              
              cont=[a*t+x_a,b*t+y_a,c*t+z_a] # Coordonnees du point d'intersection
              if Is_In(contenu,cont): # S'il y est deja, alors on a deja ajoute ce point donc on ne fait rien (on tombe deux fois sur le meme point au niveau de l'arete entre deux cellules selectionnees)
                pass
              else:
                contenu.append(cont) # On ajoute le point au contenu
                cont_ori=[cont[i]-Pt_plan[i] for i in range(0,3)]
                Pt=ProdMat(P,cont_ori)
                fichier.write("Point("+str(NB)+")={"+str(Pt[0])+", "+str(Pt[1])+", "+str(Pt[2])+", 1};"+"\n")                                                                        
                fichier2.write(str(Pt[0])+" "+str(Pt[1])+" "+str(Pt[2])+"\n")                                                                                                           
                NB+=1
        else:                                    # Idem                                                                                                                                      
          if list_test[k]*list_test[0]<0:
            x_a,y_a,z_a=Table3[Nb_precedent+4*k+1], Table3[Nb_precedent+4*k+2], Table3[Nb_precedent+4*k+3]
            x_b,y_b,z_b=Table3[Nb_precedent+1], Table3[Nb_precedent+2], Table3[Nb_precedent+3]
            norme=(x_a-x_b)**2+(y_a-y_b)**2+(z_a-z_b)**2
            if norme in length:
              a,b,c=x_b-x_a, y_b-y_a, z_b-z_a
              t=-(vect_normal[0]*x_a+vect_normal[1]*y_a+vect_normal[2]*z_a-(vect_normal[0]*Pt_plan[0]+vect_normal[1]*Pt_plan[1]+vect_normal[2]*Pt_plan[2]))/(a*vect_normal[0]+b*vect_normal[1]+c*vect_normal[2])
              cont=[a*t+x_a,b*t+y_a,c*t+z_a]
              if Is_In(contenu,cont):
                pass
              else:
                contenu.append(cont)
                cont_ori=[cont[i]-Pt_plan[i] for i in range(0,3)]
                Pt=ProdMat(P,cont_ori)
                fichier.write("Point("+str(NB)+")={"+str(Pt[0])+", "+str(Pt[1])+", "+str(Pt[2])+", 1};"+"\n")
                fichier2.write(str(Pt[0])+" "+str(Pt[1])+" "+str(Pt[2])+"\n") 
                NB+=1        
  else:
    pass
  Nb_precedent+=4*(j+1)
  num_face+=1
    
fichier2.close()
fichier3.close()
Table4=np.fromfile("Point.txt",count=-1,sep=" ")

# Ajout des points pour fermer le contour dans la direction d'extrusion
  # Obtention des coordonnees des points aux extremites droite et gauche se situant en bas de la geometrie
   # Cette acquisition n'est pas infaillible et peut etre modifie a la main a l'endroit indique
Vect=np.zeros(3)

MinMax=[]
Minl=[]
Minr=[]
min_newl=GREAT
min_newr=GREAT
alphal=GREAT
alphar=0

for i in range(0,NB):
  MinMax.append(Table4[3*i])
Delta=max(MinMax)-min(MinMax)

for i in range(0,NB):
  min_oldr=min_newr
  min_oldl=min_newl
  Vect[0],Vect[1]=Table4[3*i],Table4[3*i+1]
  Prod_scal=Vect[0]
  Prod_scal2=Vect[1]
  if Prod_scal<min(MinMax)+0.5*Delta:
    Minl.append(Prod_scal2)
    min_newl=min(Minl)
    if Prod_scal2==min_newl:
      if min_oldl==min_newl:
        if Prod_scal<alphal:
          alphal=Prod_scal
          Pt_minl=i
      else:
        alphal=Prod_scal
        Pt_minl=i
  else:
    Minr.append(Prod_scal2)
    min_newr=min(Minr)
    if Prod_scal2==min_newr:
      if min_oldr==min_newr:
        if Prod_scal>alphar:
          alphar=Prod_scal
          Pt_minr=i
      else:     
        alphar=Prod_scal
        Pt_minr=i

X_min,Y_min,Z_min=Table4[3*Pt_minl]-0.4*Delta,Table4[3*Pt_minl+1],Table4[3*Pt_minl+2] # Ajout de deux points permettant de fermer la geometrie 2D par des lignes droites
X_max,Y_max,Z_max=Table4[3*Pt_minr]+0.4*Delta,Table4[3*Pt_minr+1],Table4[3*Pt_minr+2] # Ces points peuvent directement etre modifies a la main si l'utilisateur n'est pas satisfait

  # Ajout des points pour le maillage sur les cotes et le haut
fichier.write("Point("+str(NB)+")={"+str(X_min)+", "+str(Y_min)+", "+str(Z_min)+",1};"+"\n")
fichier.write("Point("+str(NB+1)+")={"+str(X_max)+", "+str(Y_max)+", "+str(Z_max)+",1};"+"\n")
fichier.write("Point("+str(NB+2)+")={"+str(X_min)+", "+str(max(Y_min,Y_max)+H)+", "+str(Z_min)+",1};"+"\n")
fichier.write("Point("+str(NB+3)+")={"+str(X_max)+", "+str(max(Y_min,Y_max)+H)+", "+str(Z_max)+",1};"+"\n")

fichier2=open(Gmsh_points,'a')
fichier2.write(str(X_min)+" "+str(Y_min)+" "+str(Z_min)+"\n")
fichier2.write(str(X_max)+" "+str(Y_max)+" "+str(Z_max)+"\n")
fichier2.write(str(X_min)+" "+str(max(Y_min,Y_max)+H)+" "+str(Z_min)+"\n")
fichier2.write(str(X_max)+" "+str(max(Y_min,Y_max)+H)+" "+str(Z_max)+"\n")
fichier2.close()
Table4=np.fromfile("Point.txt",count=-1,sep=" ")
os.remove("Point.txt")

# Boucle for qui relie les points les plus proches les uns des autres 
  # Ne fonctionne pas toujours : 1 - Modifier le plan de coupe et le point du plan peut resoudre les problemes
  #                              2 - Modifier le fichier Maillage.geo dans un editeur de texte si le probleme est leger
condition=True
minimum=NB
list_=[i for i in range(0,NB)]
kk=0
while condition==True:
  Distance=[]
  min_new=GREAT
  i=minimum
  for j in range(0,NB):
    if j in list_:
      if i==NB:
        Distance.append((X_min-Table4[3*j])**2+(Y_min-Table4[3*j+1])**2+(Z_min-Table4[3*j+2])**2)
      else:
        Distance.append((Table4[3*i]-Table4[3*j])**2+(Table4[3*i+1]-Table4[3*j+1])**2+(Table4[3*i+2]-Table4[3*j+2])**2)
    else:
      Distance.append(GREAT)
  minimum=Distance.index(min(Distance))
  list_.remove(minimum)
  if i==NB:
    fichier.write("Line("+str(0)+")={"+str(i)+", "+str(minimum)+"};"+"\n")
  else:
    kk+=1
    fichier.write("Line("+str(kk)+")={"+str(i)+", "+str(minimum)+"};"+"\n")
  if len(list_)==0:
    break

# Ajout des lignes du haut et des cotes
fichier.write("Line("+str(kk+1)+")={"+str(minimum)+", "+str(NB+1)+"};"+"\n")
fichier.write("Line("+str(kk+2)+")={"+str(NB+1)+", "+str(NB+3)+"};"+"\n")
fichier.write("Line("+str(kk+3)+")={"+str(NB+3)+", "+str(NB+2)+"};"+"\n")
fichier.write("Line("+str(kk+4)+")={"+str(NB+2)+", "+str(NB)+"};"+"\n")


# Creation d une surface dans le fichier .geo qui sera soumis a GMSH
fichier.write("Line Loop(1) = {")
for i in range(0,kk+5):
  if i==kk+4:
    fichier.write(str(i)+"};"+"\n"+"Plane Surface(1) = {1};"+"\n")
  else:
    fichier.write(str(i)+", ")

# Ajout de noms aux limites
fichier.write('Physical Line("Boundary") = {'+str(0)+", "+str(kk+1)+", "+str(kk+2)+", "+str(kk+3)+", "+str(kk+4)+"};"+"\n")
fichier.write('Physical Line("Gas") = {')
for i in range(1,kk+1):
  if i==kk:
    fichier.write(str(i)+"};"+"\n")
  else:
    fichier.write(str(i)+", ")

fichier.write('Physical Surface("Stone") = {1};'+"\n")
fichier.write("Transfinite Line {")
for i in range(1,kk+1):
  if i==kk:
    fichier.write(str(i)+"}= "+str(Discretisation_Surf_OpenFOAM)+" Using Progression 1;"+"\n")
  else:
    fichier.write(str(i)+", ")
fichier.write("Transfinite Line {"+str(0)+"}= "+str(Discretisation_sides)+" Using Progression 1;"+"\n")
fichier.write("Transfinite Line {"+str(kk+1)+"}= "+str(Discretisation_sides)+" Using Progression 1;"+"\n")
fichier.write("Transfinite Line {"+str(kk+2)+"}= "+str(Discretisation_sides)+" Using Progression 1;"+"\n")
fichier.write("Transfinite Line {"+str(kk+4)+"}= "+str(Discretisation_sides)+" Using Progression 1;"+"\n")
fichier.write("Transfinite Line {"+str(kk+3)+"}= "+str(Discretisation_haut)+" Using Progression 1;"+"\n")
fichier.close()  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                   Correction du maillage cartesien                           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Si le maillag est cartesien, une etape de lissage est necessaire
# La procedure est decrite dans la these associee (voir le fichier "Lisez-moi.txt")
Gmsh2='CASTEM/Maillage.geo'
fichier2=open(Gmsh2,'w')
if MaillageCartesien==True:
  fichier=open(Gmsh,'r')
  Nb_line=0
  for line in fichier.readlines():
    Nb_line+=1
    if "Line" in line:                                                                       
      FirstLine=Nb_line-1
      break
  l=0
  First_Point=[]
  Second_Point=[]
  newstr=linecache.getline(Gmsh,Nb_line)
  while newstr[len(newstr)-4-l]!=" ":
    Second_Point+=newstr[len(newstr)-4-l]
    l+=1
  while newstr[len(newstr)-6-l]!="{":
    First_Point+=newstr[len(newstr)-6-l]
    l+=1
  First_Point=int("".join(reversed(First_Point)))
  Second_Point=int("".join(reversed(Second_Point)))  
  Pt1X,Pt1Y=Table4[3*First_Point],Table4[3*First_Point+1]
  Pt2X,Pt2Y=Table4[3*Second_Point],Table4[3*Second_Point+1]
  Str1="Point("+str(First_Point)+")={"+str(Pt1X)+", "+str(Pt1Y)+", 0.0, 1};\n"
  Str2="Point("+str(Second_Point)+")={"+str(Pt2X)+", "+str(Pt2Y)+", 0.0, 1};\n"
  STR=[x for x in range(0,3*kk+5)]
  STR[0]=Str1
  #STR[1]=Str2
  n=kk+5
  nbb=1
  for i in range(1,kk+5):
    l=0
    newstr=linecache.getline(Gmsh,Nb_line+i)
    Third_Point=[]
    while newstr[len(newstr)-4-l]!=" ":
      Third_Point+=newstr[len(newstr)-4-l]
      l+=1
    Third_Point=int("".join(reversed(Third_Point)))
    Pt3X,Pt3Y=Table4[3*Third_Point],Table4[3*Third_Point+1]
    Angle=((Pt2X-Pt1X)*(Pt2X-Pt3X)+(Pt2Y-Pt1Y)*(Pt2Y-Pt3Y))/(sqrt((Pt2X-Pt1X)**2+(Pt2Y-Pt1Y)**2)*sqrt((Pt3X-Pt2X)**2+(Pt3Y-Pt2Y)**2))
    if Angle>-0.2 and Second_Point<NB+1: # Alors il y a un angle droit
      Pt4X,Pt4Y=(Pt1X+Pt2X)/2.0,(Pt1Y+Pt2Y)/2.0
      Pt5X,Pt5Y=(Pt3X+Pt2X)/2.0,(Pt3Y+Pt2Y)/2.0
      Pt2X,Pt2Y=(Pt1X+4*Pt2X+Pt3X)/6.0,(Pt1Y+4*Pt2Y+Pt3Y)/6.0
      STR[nbb]="Point("+str(n)+")={"+str(Pt4X)+", "+str(Pt4Y)+", 0.0, 1};\n"
      STR[nbb+1]="Point("+str(Second_Point)+")={"+str(Pt2X)+", "+str(Pt2Y)+", 0.0, 1};\n"
      STR[nbb+2]="Point("+str(n+1)+")={"+str(Pt5X)+", "+str(Pt5Y)+", 0.0, 1};\n"
      n+=2
      nbb+=3
      Pt1X,Pt1Y=Pt5X,Pt5Y
    else:
      Pt1X,Pt1Y=Pt2X,Pt2Y
      STR[nbb]="Point("+str(Second_Point)+")={"+str(Pt2X)+", "+str(Pt2Y)+", 0.0, 1};\n"
      nbb+=1
    Pt2X,Pt2Y=Pt3X,Pt3Y
    Second_Point=Third_Point
  for i in range(0,nbb):
    fichier2.write(STR[i])
  for m in range(0,nbb-1):
    j=0
    str_ephem=""
    str_ephem2=""
    while STR[m][6+j]!=")":
      str_ephem+=str(STR[m][6+j])
      j+=1
    j=0
    while STR[m+1][6+j]!=")":
      str_ephem2+=str(STR[m+1][6+j])
      j+=1
    fichier2.write("Line("+str(m)+")={"+str(str_ephem)+", "+str(str_ephem2)+"};\n")
    
  fichier2.write("Line("+str(m+1)+")={"+str(str_ephem2)+", "+str(NB)+"};\n")
  fichier2.write("Line Loop(1) = {")
  for i in range(0,m+2):
    if i==m+1:
      fichier2.write(str(i)+"};"+"\n"+"Plane Surface(1) = {1};"+"\n")
    else:
      fichier2.write(str(i)+", ")
  # Ajout de noms aux limites
  fichier2.write('Physical Line("Boundary") = {'+str(0)+", "+str(m-2)+", "+str(m-1)+", "+str(m)+", "+str(m+1)+"};"+"\n")
  fichier2.write('Physical Line("Gas") = {')
  for i in range(1,m-2):
    if i==m-3:
      fichier2.write(str(i)+"};"+'\nPhysical Surface("Stone") = {1};\n')
    elif i<m-3:
      fichier2.write(str(i)+", ")
  fichier2.write("Transfinite Line {"+str(0)+"}= "+str(Discretisation_sides)+" Using Progression 1;"+"\n")
  fichier2.write("Transfinite Line {"+str(m-2)+"}= "+str(Discretisation_sides)+" Using Progression 1;"+"\n")
  fichier2.write("Transfinite Line {"+str(m-1)+"}= "+str(Discretisation_sides)+" Using Progression 1;"+"\n")
  fichier2.write("Transfinite Line {"+str(m)+"}= "+str(Discretisation_haut)+" Using Progression 1;"+"\n")
  fichier2.write("Transfinite Line {"+str(m+1)+"}= "+str(Discretisation_sides)+" Using Progression 1;"+"\nTransfinite Line {") 
  for i in range(1,m-2):
    if i==m-3:
      fichier2.write(str(i)+"}="+str(Discretisation_Surf_OpenFOAM)+" Using Progression 1;"+"\n")
    elif i<m-3:
      fichier2.write(str(i)+", ")   
  fichier.close()
  fichier2.close() 
else:
  os.rename('CASTEM/Maillages.geo','CASTEM/Maillage.geo') 
  
  
