from string import *
import os
import sys

def aout(bif1):
    for s in bif1:
	#on utilise bash pour faire un fichier de parametres

	a= "sed -e \"s/bif1/"+replace(s,'\n','')+"/\""+" params>tata"
	os.system(a)

        cmd= "./a.out<"+"tata"
	handle = os.popen(cmd,'r',1)
	for line in handle:
	        print line,
	handle.close()

        renome= "mv traj.dat traj."+s
	os.system(renome)
	## os.system("cat peak1.dat >> peak1.bif")
	## os.system("cat peak2.dat >> peak2.bif")
	## os.system("cat seuil.dat >> seuil.bif")

os.chdir('hb_ri')
print os.getcwd()
bif1=open("bif1")
aout(bif1)
bif1.close()

os.chdir('../hb_rs')
print os.getcwd()
bif1=open("bif1")
aout(bif1)
bif1.close()

os.chdir('../sb_ri')
print os.getcwd()
bif1=open("bif1")
aout(bif1)
bif1.close()

os.chdir('../sb_rs')
print os.getcwd()
bif1=open("bif1")
aout(bif1)
bif1.close()

os.chdir('../gupta')
print os.getcwd()
bif1=open("bif1")
aout(bif1)
bif1.close()

