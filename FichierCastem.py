import numpy as np
from math import *
import os
import linecache
import sys
import csv
from Mailleur import Name_surf,Version_castem,vect_normal,Pt_plan,Pas,run,rupelian
GREAT=1e10
# Les caracteristiques du calcul Castem sont a modifier a la main

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                           Temperature initiale                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Choisi la temperature initiale du gaz par defaut et non la temperature initiale des surfaces (autrement bugs possibles avec $internalField)
PATH_ini="0/T" 
Temp0=open(PATH_ini, 'r')
inside=Temp0.read()
N_inside=inside.find('internalField')
inside_2=inside[N_inside:]
N_inside2=inside_2.find(';')
nt=0
nbt=""
tar=""
while tar!=" ":
  nt+=1
  tar=inside_2[N_inside2-nt]
  nbt+=tar
  
Temp_ini=float("".join(reversed(nbt)))-273.15


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  Obtention du numero de la premiere ligne dans les fichiers de temperature   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
PasDt=[]
Temp_line=[] # Liste des lignes de la premiere temperature pour chaque fichier X/T avec X le temps de calcul enregistre dans OpenFOAM
Garbage=[]


for file in os.listdir('.'):
  try:
    file=float(file)
    if file!=0:
      if file==int(file):
        Garbage.append(int(file))
      else:
        Garbage.append(file)
  except ValueError:
    pass

for i in range(0,len(Garbage)):
  PasDt.append(min(Garbage)) # Contient tous les chemins vers les fichiers X avec X le temps de calcul enregistre dans OpenFOAM
  Garbage.remove(min(Garbage))

PATH=[]
q=0

for dt in PasDt:
  Temp_line_surf=np.zeros(len(Name_surf)) # Pour un dossier X/T, contient la premiere ligne de la premiere temperature pour chaque surface
  PATH.append(str(dt)+"/T")
  for i in range(0,len(Name_surf)):
    nligne=0
    Temp=open(PATH[q],'r')   # Ouverture du fichier texte X/T
    for ligne in Temp:
      nligne+=1
      if Name_surf[i] in ligne:
        Surf_ligne=nligne
        break
    Temp.close()
    Temp=open(PATH[q],'r')
    lines=Temp.readlines()[Surf_ligne:]
    nline=0
    for line in lines:
      nline+=1
      if "value " in line:
        Temp_line_surf[i]=nline+2+Surf_ligne
        break
    Temp.close()
  Temp_line.append(Temp_line_surf)
  #Temp.close()
  q+=1


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  Ecriture du debut du fichier dgibi et recuperation du maillage dans Cast3m  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
Castem='CASTEM/Thermo_mecanique.dgibi'
fichier=open(Castem,'w')
# Definition du type d'element utilise pour le calcul
fichier.write("OPTI DIME 2 ELEM TRI3; \n")
# Contraintes planes
fichier.write("*Hypotheses des deformations planes\nOPTI 'MODE' 'PLAN' 'DEFO';\n") 
# Lecture du maillage
MaillageUNV="CASTEM/Maillage.unv"
fichier.write("TAB1=LIRE 'UNV' 'Maillage.unv';\n")
# Recuperation de la surface
fichier.write("SU=TAB1.'Stone';\n")
# Recuperation des elements a la limite
fichier.write("GAS=TAB1.'Gas';"+"\n")
fichier.write("SIDES=TAB1.'Boundary';"+"\n")
# Nombre d'elements contenant un bord confondu avec la limite Gas
Nb_line=0
with open(MaillageUNV,"r") as f:
  for line in f.readlines():
    Nb_line+=1
    if "Gas" in line:                                                                             
      Nb_Gas=Nb_line-1
      break
chaine=linecache.getline(MaillageUNV,Nb_Gas)
n=0
NB_Gas=""
while chaine[len(chaine)-2-n]!=" ":
  NB_Gas+=chaine[len(chaine)-2-n]
  n+=1
NB_Gas=int("".join(reversed(NB_Gas)))

# Nombre d'elements contenant un bord confondu avec la limite Boundary
Nb_line=0
with open(MaillageUNV,"r") as f:
  for line in f.readlines():
    Nb_line+=1
    if "Boundary" in line:                                                                       
      Nb_Bound=Nb_line-1
      break
chaine2=linecache.getline(MaillageUNV,Nb_Bound)
n=0
NB_Bound=""
while chaine2[len(chaine2)-2-n]!=" ":
  NB_Bound+=chaine2[len(chaine2)-2-n]
  n+=1
NB_Bound=int("".join(reversed(NB_Bound)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Recuperation des coordonnees des lignes du maillage et calcul de leur centre #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Copie de Thermo_mecanique.dgibi
fich=open("Poubelle.dgibi","w")
fich.write("OPTI DIME 2 ELEM TRI3;\n")
# Lecture du maillage
fich.write("TAB1=LIRE 'UNV' '"+str(MaillageUNV)+"';\n")
# Recuperation des elements a la limite
fich.write("GAS=TAB1.'Gas';\n")                                                                          

# Recuperation des barycentres de chaque ligne dans un fichier DETRITUS.csv
fich.write("NB_GAS=NBEL GAS;\n")
fich.write("REPE I NB_GAS; \n L=GAS ELEM &I; \n P=BARY L; \n R1 R2=COOR P; \n COND = &I <EG 1.; \n") 
fich.write(" SI COND; \n  LIST1=PROG R1 R2; \n SINON; \n  LIST2=PROG R1 R2; \n  LIST1=LIST1 ET LIST2; \n FINSI; \nFIN I; \n")
fich.write("OPTI 'SORT' DETRITUS; \n SORT 'EXCE' LIST1 'SEPA' 'TABU'; \nFIN;")
fich.close()
executable="castem"+str(Version_castem)+" Poubelle.dgibi"
os.system(executable) 
os.remove("Poubelle.dgibi")
os.remove("Poubelle.trace")
Barycentre=np.zeros((NB_Gas+1,3)) # Stockage des barycentres dans ce vecteur dans l'ordre des lignes
ncount=0

with open('DETRITUS.csv', 'rb') as f:
  reader=csv.reader(f,delimiter='\t')
  for row in reader:
    if ncount!=0:
      if ncount % 2==1:
        Barycentre[int((ncount-1)/2.0),0]=row[0]
      elif ncount % 2==0:
        Barycentre[int((ncount-1)/2.0),1]=row[0]
    ncount+=1
os.remove("DETRITUS.csv")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  Ecriture de la liste des temps calcules par OpenFOAM dont on a les donnees  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fichier.write("LTime_0=PROG 0 ")
lenstr="LTime_0=PROG 0 "
reset=0
l=0
for npath in range(0,len(PATH)):
  fichier.write(str(PasDt[npath])+" ")
  lenstr+=str(PasDt[npath])+" "
  if npath==len(PATH)-1:
    fichier.write(";\n")
    break
  if len(lenstr)>63:
    fichier.write("\n")
    lenstr=""
    reset+=1
  if reset==6:
    l+=1
    fichier.write(";\n"+"LTime_"+str(l)+"=PROG ")
    reset=0
    lenstr="LTime_"+str(l)+"=PROG "

fichier.write("LTime= ")
strl="LTime= "
i=0
for n in range(0,l+1):
  if n<l:
    strl+="LTime_"+str(n)+" ET "
    fichier.write("LTime_"+str(n)+" ET ")
  else:
    fichier.write("LTime_"+str(n)+";\n")
  if len(strl)>61:
    fichier.write("\n")
    strl=""


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#             Recuperation des barycentres des faces d'OpenFOAM                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
Table=np.fromfile("Castem.txt",count=-1,sep=" ")
os.remove("Castem.txt")
ListFaceIntersect=[] # Contient toutes les faces intersectees par le plan de coupe
ListBarycentre=np.zeros((len(Table)/5,3)) # Contient le barycentre correspondant a la face de ListFaceIntersect
ListFaceSurface=[] # Contient le numero de la surface auquel appartient la face correspondante dans ListFaceIntersect
for i in range(0,len(Table)/5):
  ListFaceSurface.append(int(Table[5*i]))
  ListFaceIntersect.append(int(Table[5*i+1]))
  ListBarycentre[i,0],ListBarycentre[i,1],ListBarycentre[i,2]=Table[5*i+2],Table[5*i+3],Table[5*i+4]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   Listes de temperatures en fonction du temps pour chaque face d'OpenFOAM    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
number=0
for ii in ListFaceIntersect:
  fichier.write("LT"+str(number)+"_0=PROG "+str(Temp_ini)+" ")
  lenstr="LT"+str(number)+"_0=PROG "+str(Temp_ini)+" "
  reset=0
  l=0
  for j in range(0,len(PATH)):
    ficTemp=open(PATH[j],'r')
    newstr=str(ficTemp.readlines()[int(Temp_line[j][ListFaceSurface[number]-1])+ii])
    newstr=str(float(newstr.rstrip())-273.15)  
    fichier.write(newstr+" ")
    lenstr+=newstr+" "
    ficTemp.close()
    if j==len(PATH)-1:
      fichier.write(";\n")
      break
    if len(lenstr)>63:
      fichier.write("\n")
      lenstr=""
      reset+=1
    if reset==6:
      l+=1
      fichier.write(";\n"+"LT"+str(number)+"_"+str(l)+"=PROG ")
      reset=0
      lenstr="LT"+str(number)+"_"+str(l)+"=PROG "
  fichier.write("LT"+str(number)+"= ")
  strl="LT"+str(number)+"= "
  for n in range(0,l+1):
    if n<l:
      strl+="LT"+str(number)+"_"+str(n)+" ET "
      fichier.write("LT"+str(number)+"_"+str(n)+" ET ")
    else:
      fichier.write("LT"+str(number)+"_"+str(n)+";\n")
    if len(strl)>59:
      fichier.write("\n")
      strl=""
  fichier.write("EVT"+str(number)+"=EVOL 'MANU' 'Temps' LTime 'Coef' LT"+str(number)+";\n")
  number+=1
fichier.write("\n")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Association des lignes de la limite GAS aux faces d OpenFOAM          #
#      Attribution du champ de temperature a chaque ligne du maillage 2D       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
STRGAS=[]
Veclength=[]
for k in range(0,len(Table)/5):
	STRGAS.append("LI"+str(k)+"=LECT")
	Veclength.append("LI"+str(k)+"=LECT")

# Pour chaque ligne, on choisit la temperature de la face dont le barycentre est le plus proche de cette ligne
for i in range(1,NB_Gas+1):
	X,Y,Z=Barycentre[i,0], Barycentre[i,1], 0
	Dist=[(X-ListBarycentre[k,0])**2+(Y-ListBarycentre[k,1])**2+(Z-ListBarycentre[k,2])**2 for k in range(0,len(Table)/5)]
	ind=Dist.index(min(Dist))	
	if len(Veclength[ind])<63:
		STRGAS[ind]+=" "+str(i)
		Veclength[ind]+=" "+str(i)
	else:
		Veclength[ind]=str(i)
		STRGAS[ind]+="\n "+str(i)

i=0
lol=0
for k in range(0,len(Table)/5):
	if STRGAS[k]=="LI"+str(k)+"=LECT":
		continue
	else: 
		STRGAS[k]+=";\n"
		l=STRGAS[k].count('\n')
		if l>6:
			recur=0
			while 1==1:
				j=400
				if len(STRGAS[k])<j:
					AnotherChaine="LI"+str(len(Table)/5+lol)+"=LECT "+str(STRGAS[k])+";\n"		
					fichier.write(AnotherChaine)
					fichier.write("L"+str(i)+" = GAS ELEM LI"+str(len(Table)/5+lol)+";\n")
					fichier.write("B0"+str(i)+" = BLOQ L"+str(i)+" 'T';\n")
  					fichier.write("D0"+str(i)+" = DEPI B0"+str(i)+" 1.;\n")
					fichier.write("C0"+str(i)+" = CHAR 'TIMP' D0"+str(i)+" EVT"+str(k)+";\n")					
					i+=1
					lol+=1
					break
				else:
					while STRGAS[k][j]!='\n':
						j+=1
					j-=50
					while STRGAS[k][j]!=" ":
						j+=1
					if recur==0:
						AnotherChaine=STRGAS[k][0:j]+";\n"
						fichier.write(AnotherChaine)
						fichier.write("L"+str(i)+" = GAS ELEM LI"+str(k)+";\n")
					else:
						AnotherChaine="LI"+str(len(Table)/5+lol)+"=LECT "+STRGAS[k][0:j]+";\n"		
						fichier.write(AnotherChaine)
						fichier.write("L"+str(i)+" = GAS ELEM LI"+str(len(Table)/5+lol)+";\n")
						lol+=1
					fichier.write("B0"+str(i)+" = BLOQ L"+str(i)+" 'T';\n")
  					fichier.write("D0"+str(i)+" = DEPI B0"+str(i)+" 1.;\n")
					fichier.write("C0"+str(i)+" = CHAR 'TIMP' D0"+str(i)+" EVT"+str(k)+";\n")
					i+=1
					recur+=1
					STRGAS[k]=STRGAS[k].replace(STRGAS[k][0:j],'')
		else:
			fichier.write(STRGAS[k])
			fichier.write("L"+str(i)+" = GAS ELEM LI"+str(k)+";\n")
			fichier.write("B0"+str(i)+" = BLOQ L"+str(i)+" 'T';\n")
			fichier.write("D0"+str(i)+" = DEPI B0"+str(i)+" 1.;\n")
			fichier.write("C0"+str(i)+" = CHAR 'TIMP' D0"+str(i)+" EVT"+str(k)+";\n")
			i+=1

II=i
# Boucle attribuant un nom a chaque ligne limitrophe
for i in range(1,NB_Bound+1):
  fichier.write("LB_"+str(i)+" = SIDES ELEM "+str(i)+";"+"\n")

# Puisque les versions anciennes de Castem n'acceptent pas plus de 72 caracteres par ligne, il est necessaire d'utiliser plusieurs lignes pour la directive
# De plus, chaque directive ne doit pas exceder 500 caracteres. Il est donc necessaire d'utiliser plusieurs directives pour creer chaque liste
# C'est l'objet de toute la suite jusqu'au prochain bloc
def fwrite3(i,k,f,maxx):
	str1="B"+str(f) + str(i) +"="
	m=k
	while len(str1)<61 and m<=maxx:
		str1+="B"+str(f-1)+str(m)+" ET "
		m=m+1
	res=str1[:-4]
	res=res+";"+"\n"
	return res
	
f=0
nb=2
y=II
while nb>1:
	nb=0
	g=0
	o=0
	f=f+1
	maxx=y-1
	while o<=y-1:
		nb=nb+1
		sl=fwrite3(g,o,f,maxx)
		fichier.write(str(sl))
		g=g+1
		z=3
		cd=""
		while sl[len(sl)-z]!="B":
			cd+=str(sl[len(sl)-z])
			z=z+1
		cdd=cd[::-1]
		if f<10:
			o=int(cdd[1:])+1
		elif f<100 and f>9:
			o=int(cdd[2:])+1
	y=g

def fwrite4(i,k,f,maxx):
	str1="C"+str(f) + str(i) +"="
	m=k
	while len(str1)<61 and m<=maxx:
		str1=str1+"C"+str(f-1)+str(m)+" ET "
		m=m+1
	res=str1[:-4]
	res=res+";"+"\n"
	return res
	
f=0
nb=2
y=II
while nb>1:
	nb=0
	g=0
	o=0
	f=f+1
	maxx=y-1
	while o<=y-1:
		nb=nb+1
		sl=fwrite4(g,o,f,maxx)
		fichier.write(str(sl))
		g=g+1
		z=3
		cd=""
		while sl[len(sl)-z]!="C":
			cd=cd+str(sl[len(sl)-z])
			z=z+1
		cdd=cd[::-1]
		if f<10:
			o=int(cdd[1:])+1
		elif f<100 and f>9:
			o=int(cdd[2:])+1
	y=g
	
# Conditions limites en deplacement
fichier.write("BLMX=BLOQ SIDES 'UX';\nBLMY=BLOQ SIDES 'UY';\n")
# Sauvegarde des elements et des conditions limites dans Cast3m
fichier.write("\n*NOM DU FICHIER DE SAUVEGARDE\nOPTI 'SAUV' 'ThermoMeca.sauv';\n*ECRITURE DES FICHIERS\nSAUV;\n")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                  Ecriture des caracteristiques materiau                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# A modifier selon le materiau
if rupelian==True:
	fichier.write("\n* ========================================================== *\n*												Materiau 			   									  *\n* ========================================================== *\n*Caracteristique du materiau"+"\n"+"RHOMAT=EVOL MANU 'T'(PROG 20. 1300.)"+"\n"+"'RHO' (PROG 1675. 1675.);"+"\n"+"CAPAMAT=EVOL MANU 'T' (PROG 20. 95. 100. 110. 115. 300. 400. 1300.)"+"\n"+"'C' (PROG 669. 669. 10000. 10000. 669. 669. 669. 669.);"+"\n"+"CONDUMAT=EVOL MANU 'T' (PROG 20. 100. 200. 300. 400. 500. 1200.) "+"\n"+"'K' (PROG 0.74 0.54 0.43 0.37 0.31 0.25 0.25);"+"\n"+"ALPHAMAT=EVOL MANU 'T' (PROG 20. 300. 1200.) "+"\n"+"'ALPH' (PROG 2.e-6 1.6E-5 1.6E-5);"+"\n"+"YOUNGMAT=EVOL MANU 'T' (PROG 20. 250. 600. 1200.) "+"\n"+"'YOUN' (PROG 2.5E9 2.0E9 2.0E9 1.0E6);"+"\n"+"NUMAT=0.3;"+"\n"+"SIGYMAT=150.E6;"+"\n")
else:
	fichier.write("\n* ========================================================== *\n*												Materiau 			   									  *\n* ========================================================== *\n*Caracteristique du materiau"+"\n"+"RHOMAT=EVOL MANU 'T'(PROG 20. 1300.)"+"\n"+"'RHO' (PROG 2180. 2180.);"+"\n"+"CAPAMAT=EVOL MANU 'T' (PROG 20. 95. 100. 110. 115. 300. 400. 1300.)"+"\n"+"'C' (PROG 894. 894. 10000. 10000. 894. 894. 894. 894.);"+"\n"+"CONDUMAT=EVOL MANU 'T' (PROG 20. 100. 200. 300. 400. 500. 1200.) "+"\n"+"'K' (PROG 1.83 1.34 1.07 0.92 0.76 0.61 0.61);"+"\n"+"ALPHAMAT=EVOL MANU 'T' (PROG 20. 250. 1200.) "+"\n"+"'ALPH' (PROG 4.e-6 1.0E-5 1.0E-5);"+"\n"+"YOUNGMAT=EVOL MANU 'T' (PROG 20. 250. 600. 1200.) "+"\n"+"'YOUN' (PROG 2.5E10 9.0E9 9.0E9 1.0E6);"+"\n"+"NUMAT=0.3;"+"\n"+"SIGYMAT=150.E6;"+"\n")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#            Procedure PASAPAS pour la resolution thermo-mecanique             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fichier.write("\n* ========================================================== *\n*									Modelisation du probleme	      				   *\n* ========================================================== *\n*MODELE THERMIQUE (CONDUCTION) MATERIAU UNIFORME ET CONSTANT\n*MODELE MECANIQUE MATERIAU ELASTIQUE ISOTROPE\n\nMOT=MODE SU 'THERMIQUE' 'ISOTROPE' 'CONDUCTION';\nMAT=MATE MOT 'K' CONDUMAT 'C' CAPAMAT 'RHO' RHOMAT;\nMOM=MODE SU MECANIQUE ELASTIQUE ISOTROPE;\nMAM=MATE MOM 'YOUN' YOUNGMAT 'NU' NUMAT 'ALPH' ALPHAMAT;\n")
fichier.write("\n* ========================================================== *\n* Resolution de la thermo-mecanique par la procedure PASAPAS *\n* ========================================================== *\n*Champ de temperatures initiales\nT_INI=MANU 'CHPO' SU 1 'T' "+str(Temp_ini)+";\n")

fichier.write("* Instant final du calcul\nTPSFIN="+str(PasDt[len(PasDt)-1])+";\nPASCAL="+str(Pas)+";\nPASAUV=10*PASCAL;\n")

fichier.write("*Definition de la table d'arguments\n*Procedure PASAPAS\nTAB2=TABL;\nTAB2. 'MODELE'= MOT ET MOM;\nTAB2. 'CARACTERISTIQUES' = MAT ET MAM;\nTAB2 . 'CHARGEMENT'=C"+str(f)+str(0)+";\nTAB2. 'BLOCAGES_THERMIQUES' = B"+str(f)+str(0)+";\nTAB2.'BLOCAGES_MECANIQUES'=BLMX ET BLMY;\nTAB2 . 'TEMPS_CALCULES' = PROG 0. 'PAS' PASCAL TPSFIN ;\n*TAB2 . 'AUTOMATIQUE' = VRAI;\n*TAB2.'AUTOCRIT'=0.001;\n*TAB2.'AUTOPAS'=100;\n*TAB2.'AUTORESU'=5;\n*TAB2 . 'PAS_AJUSTE' = VRAI;\nTAB2. 'TEMPS_SAUVES'= PROG 0. 'PAS' PASAUV TPSFIN;\nTAB2 . 'TEMPERATURES' = TABL ;\nTAB2 . 'TEMPERATURES' . 0 = T_INI ;\nTAB2 . 'CELSIUS' =VRAI;\n\n*APPEL PASAPAS\nPASAPAS TAB2 ;\n")

fichier.write("\n*NOM DU FICHIER DE SAUVEGARDE\nOPTI 'SAUV' 'ThermoMeca.sauv';\n*ECRITURE DES FICHIERS\nSAUV;\n")
fichier.write("OPTI 'REST' 'ThermoMeca.sauv';\nREST;\n\n")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#             Ajout du trace du champ de temperature a chaque pas              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fichier.write("N1=DIME(TAB2.'TEMPERATURES');\nREPE B1 N1;\nT_I=TAB2.'TEMPERATURES' .(&B1-1);\nTPS_I=(&B1-1)*PASCAL;\nMOT_I=CHAI '[4] Temperatures au temps' TPS_I;\n*TRAC T_I SURF1 'TITR' MOT_I (PROG 25. 'PAS' 2.5 300.);\nFIN B1;\nTRAC T_I SU 'TITR' 'Temperatures au temps final';\n\n")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#              Ajout du trace de la deformee a la fin du calcul                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fichier.write("N2=DIME (TAB2.'DEPLACEMENTS');\nU=TAB2.'DEPLACEMENTS'.(N2-1);\nDEF=DEFO SU U 150. 'BLEU';\nDEF_INI=DEFO SU U 0.;\nTRAC (DEF_INI ET DEF) 'TITR' 'Deformee';\n\n")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#            Ajout du champ de contrainte a chaque pas enregistre              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fichier.write("*PROCEDURE QUI TRACE N'IMPORTE QUEL CHAMP A TOUS LES TEMPS SAUVES\nDEBP @PTT T1*'TABLE' MOT1*'MOT' MOT2/'MOT' LIS1/'LISTREEL';\nNPAS=DIME(T1.'TEMPS');\nMO=EXTR(T1.'MODELE') 'FORM' 'MECANIQUE';\nMA=REDU(T1.'CARACTERISTIQUES') MO;\nMAIL1=EXTR MO 'MAIL';\nCONT1=CONT MAIL1;\nREPE B1 NPAS;\nTPS1=T1.'TEMPS'.(&B1-1);\nDEPL1=T1.'DEPLACEMENTS'.(&B1-1);\nDEF1=DEFO MAIL1 DEPL1 150.;\nCHAMP1=T1.MOT1.(&B1-1);\nTIT1=CHAI MOT1;\nSI(NEG (TYPE MOT2) 'ANNULE');\nCHAM1=EXCO CHAM1 MOT2;\nTIT1=CHAI TIT1 ', composante ' MOT2;\nFINSI;\nTIT1=CHAI TIT1 ', au temps ' TPS1;\nSI(EGA (TYPE LIS1) 'ANNULE');\nTRAC CHAM1 MO DEF1 CONT1 15 'TITR' TIT1;\nSINON;\nTRAC CHAM1 MO DEF1 CONT1 'TITR' TIT1 LIS1;\nFINSI;\nFIN B1;\nFINP;\n")
fichier.write("*CONTRAINTES ET DEFORMATIONS PLASTIQUES CUMULEES\n*@PTT TAB2 'CONTRAINTES' (PROG 0. 'PAS' 10.E6 160.E6);\n\n")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                           Critere Rubefaction                                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fichier.write("* ========================================================== *\n*            Fonction rubefaction          *\n* ========================================================== *\n*PROCEDURE DE RUBEFACTION\nDEBP @RUBE TAB*'TABLE';\nN3=DIME(TAB.'TEMPERATURES');\nTCH=CHAN 'CHAM' T_INI MOT;\nTCH3=TCH MASQ SUPERIEUR 2000;\nTCG3=TCH MASQ SUPERIEUR 2000;\nREPE B3 N3;\nTPS_I=(&B3-1)*PASCAL;\nT_J=TAB.'TEMPERATURES'.(&B3-1);\nTCHA=CHAN 'CHAM' T_J MOT;\nTCH1=TCHA/250.;\nTCG1=TCHA/350.;\nTCH2=TCH1 MASQ SUPERIEUR 1;\nTCG2=TCG1 MASQ SUPERIEUR 1;\nTCH3=TCH3+TCH2;\nTCG3=TCG3+TCG2;\nFIN B3;\nFINP TCH3 TCG3;\nRUB_R RUB_G=@RUBE TAB2;\n")
fichier.write("RUB_R=PASCAL*RUB_R;\nRUB_G=PASCAL*RUB_G;\nRED=CHAN 'CHPO' MOT RUB_R;\nGREY=CHAN 'CHPO' MOT RUB_G;\n")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                            Procedures paraview                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fichier.write("ntot = dime TAB2.TEMPERATURES;\na = 0;\n\nREPE I ntot;\n  DD = TAB2.TEMPERATURES.a;\n  OPTI SORT  (chain'T'a'.inp');\n  SORT AVS DD SU;\n  a = a + 1;\nFIN I;\n")
fichier.write("\na = 0;\nREPE I ntot;\n  DD = TAB2.CONTRAINTES.a;\n  OPTI SORT  (chain'Sigma'a'.inp');\n  SORT AVS DD SU;\n  a = a + 1;\nFIN I;\n")
fichier.write("OPTI SORT (chain'Rubefaction');\nSORT AVS RED SU;\nOPTI SORT (chain'CouleurGrise');\nSORT AVS GREY SU;\n")


fichier.write("FIN;")
fichier.close()

os.remove("Mailleur.pyc")

if run==True:
  os.chdir("CASTEM")
  executable2="castem"+str(Version_castem)+" Thermo_mecanique.dgibi"
  os.system(executable2) 
else:
  pass

