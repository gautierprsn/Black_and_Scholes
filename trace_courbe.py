##PUT
#Crank Nicholson
f=open(r"C:\Users\gauti_000\Documents\Travail\ecole ingé\2eme annee\PAP\Projet\projet\result_put.txt","r")
s=f.readlines()

x=[0 for i in range(int(len(s)/2))]
j=0

for i in range(len(x)):
    chaine=''
    while(s[len(s)-2][j] != ";" and j<len(s[len(s)-2])):
        chaine=chaine+s[len(s)-2][j]
        j=j+1
    x[i]=float(chaine)
    j=j+1
    
#Differences finies
f=open(r"C:\Users\gauti_000\Documents\Travail\ecole ingé\2eme annee\PAP\Projet\projet\result_put_DF.txt","r")
s1=f.readlines()

x1=[0 for i in range(int(len(s1)/2)-1)]
j=0

for i in range(len(x1)):
    chaine=''
    while(s1[len(s1)-2][j] != ";" and j<len(s1[len(s1)-2])):
        chaine=chaine+s1[len(s1)-2][j]
        j=j+1
    x1[i]=float(chaine)
    j=j+1

erreur=[0 for i in range(len(x1))]
for i in range(100):
    erreur[i]=x[i]-x1[i]

import matplotlib.pyplot as plt
plt.figure()
plt.plot(x,label="Crank Nicholson")
plt.plot(x1,label="Différences Finies")
plt.title("Allure d'un put avec la methode Crank Nicholson et la methode Différences Finies")
plt.legend()
plt.show()

plt.figure()
plt.plot(erreur,label="écart entre les 2 méthodes pour un put")
plt.title("Erreur pour un put")
plt.legend()
plt.show()


##CALL
f=open(r"C:\Users\gauti_000\Documents\Travail\ecole ingé\2eme annee\PAP\Projet\projet\result_call.txt","r")
s2=f.readlines()

x2=[0 for i in range(int(len(s2)/2))]
j=0

for i in range(len(x2)):
    chaine=''
    while(s2[len(s2)-2][j] != ";" and j<len(s2[len(s2)-2])):
        chaine=chaine+s2[len(s2)-2][j]
        j=j+1
    x2[i]=float(chaine)
    j=j+1
    
    
f=open(r"C:\Users\gauti_000\Documents\Travail\ecole ingé\2eme annee\PAP\Projet\projet\result_call_DF.txt","r")
s3=f.readlines()

x3=[0 for i in range(int(len(s3)/2))]
j=0

for i in range(len(x3)):
    chaine=''
    while(s3[len(s3)-2][j] != ";" and j<len(s3[len(s3)-2])):
        chaine=chaine+s3[len(s3)-2][j]
        j=j+1
    x3[i]=float(chaine)
    j=j+1

erreur=[0 for i in range(len(x2))]
for i in range(len(erreur)):
    erreur[i]=x2[i]-x3[i]


import matplotlib.pyplot as plt
plt.figure()
plt.plot(x2,label="Crank Nicholson")
plt.plot(x3,label="Différences Finies")
plt.title("Allure d'un call avec la methode Crank Nicholson et la methode Différences Finies")
plt.legend()
plt.show()

plt.figure()
plt.plot(erreur,label="écart entre les 2 méthodes pour un call")
plt.title("Erreur pour un call")
plt.legend()
plt.show()