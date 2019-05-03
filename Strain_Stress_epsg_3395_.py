
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 14 21:16:39 2014

@author: astronomia
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import pyproj as pyproj

prjMerc=pyproj.Proj("+init=EPSG:3395")



f=open('C:/Users/JulianJar92/source/repos/Strain_Stress_epsg(3395)/entrada_datos.txt','r')
s=f.readlines()
f.close()
filesDatos=s[1].replace('\n', '').split()
print filesDatos
nombresFiles=s[3].replace('\n', '').split()
print nombresFiles

sufijo='_Dela_Excluir_LPAL_MAS1_'
excluir=['LPAL','MAS1']
#sufijo='_Dela'
#excluir=[]

dias=float(s[5].replace('\n', ''))
pasoyuni=s[7].replace('\n', '').split()
pasoametros = float(pasoyuni[0])
unidades = pasoyuni[1]
horver=s[9].replace('\n', '').split()
hori=float(horver[0])
vert=float(horver[1])
escalas=s[11].replace('\n', '').split()
hscale=float(escalas[0])
vscale=float(escalas[1])
limites=s[13].replace('\n', '').split()
xmin=float(limites[0])
xmax=float(limites[1])
ymin=float(limites[2])
ymax=float(limites[3])
escalas=s[15].replace('\n', '').split()
xescala=float(escalas[0])
yescala=float(escalas[1])
escalas=s[17].replace('\n', '').split()
strainscale=float(s[17])
escalas=s[19].replace('\n', '').split()
tramos=float(s[19])
filemapa=s[21].replace('\n', '')

for contador in range(len(filesDatos)):
    filedatos=filesDatos[contador]
    nomfil=nombresFiles[contador]+sufijo
    print('fichero de datos: ' + filedatos)
    fd = open(filedatos,'r')    
    d = np.loadtxt(fd,
                   delimiter=',',
                   dtype={ 'names':('PAIS','vertex_index','Long','Lat','Evel(mm/año)','Nvel  (mm/año)','Uvel  (mm/año)','Huso','Hemisferio','X(utm)','Y(utm)','Estacion','PAIS2'),
                          'formats':('S9','int','float','float','float','float','float','int','S9','float','float','S4','S9')},
                   comments='%')
    fd.close()


    sites=np.array([ele[0] for ele in d if (ele[0] not in excluir)]).transpose()
    print sites
    a=np.array([[ele[i] for ele in d if (ele[11] not in excluir)] for i in [2,3,4,5,6]]).transpose()
    print a

    #%velocidades desplazamiento en mm 
    for i in range(len(a)):
        a[i,0],a[i,1]=prjMerc(a[i,0],a[i,1])
        a[i,2]=a[i,2]*365*1000/pasoametros/dias
        a[i,3]=a[i,3]*365*1000/pasoametros/dias
        a[i,4]=a[i,4]*365*1000/pasoametros/dias

    puntos=np.array([[a[i,0],a[i,1]] for i in range(len(a))])
    print puntos
    tri=Delaunay(puntos)

    n=len(tri.simplices)

    fidout = open('Dela_Strain_Stress_' + nomfil + '.txt','w+')
    fidout.write('X,Y,ESTE,NORTE,EMAX,EMIN,THETA,ROTACION,DILATACION,CIZALLA,MAXGEO,EMAXABS,EMINABS,THETA90,THETA180,THETA270\n')
    X0=np.zeros([n])
    Y0=np.zeros([n])
    s=0
    D=np.empty(shape=[0,14]) 
    x=np.zeros(n)
    y=np.zeros(n)
    for i in range(n):
        x[0]=puntos[tri.simplices[i,0]][0]
        y[0]=puntos[tri.simplices[i,0]][1]    
        x[1]=puntos[tri.simplices[i,1]][0]
        y[1]=puntos[tri.simplices[i,1]][1]    
        x[2]=puntos[tri.simplices[i,2]][0]
        y[2]=puntos[tri.simplices[i,2]][1]    
        x0=(x[0]+x[1]+x[2])/3
        y0=(y[0]+y[1]+y[2])/3
    
        s=0
        A=np.zeros([6,6])
        V=np.zeros([6,1])
        s=0
        for k in [0,1,2]:
            A[(s+1)*2-2,0]=x[k]-x0
            A[(s+1)*2-2,1]=y[k]-y0
            A[(s+1)*2-2,2]=1
            A[(s+1)*2-1,3]=x[k]-x0
            A[(s+1)*2-1,4]=y[k]-y0
            A[(s+1)*2-1,5]=1
            V[(s+1)*2-2]=a[tri.simplices[i,k],2]/1000
            V[(s+1)*2-1]=a[tri.simplices[i,k],3]/1000
            s=s+1
        N=np.dot(A.transpose(),A)
        V=np.dot(A.transpose(),V)
        B=np.dot(np.asmatrix(N)**-1,V)
        exx=B[0]
        eyy=B[4]
        exy=(B[1]+B[3])/2
        rotation=(B[1]-B[3])/2
        dilatacion=exx+eyy
        cizalla1=exx-eyy
        cizalla2=2*exy
        maxcizalla=np.sqrt(cizalla1**2+cizalla2**2)
        thetastrain=.5*np.arctan(-cizalla2/cizalla1) + np.pi/2
        emax=(dilatacion+maxcizalla)/2
        emin=(dilatacion-maxcizalla)/2
        maxgeo=max(np.abs(emax),np.abs(emin))
        #este, norte, elevacion, delta gamma emax emin fi-strain max-defor-geo inclinacion fi-inclinacion
        D=np.append(D,[[x0, y0, B[2]*1000, B[5]*1000, 0, dilatacion, maxcizalla, emax, emin, thetastrain,  maxgeo, 0, 0, rotation]],axis=0)

        fidout.write('%-8.2f,%-8.2f,%-8.5f,%-8.5f,%-8.5f,%-8.5f,%-8.5f,%-8.5f,%-8.5f,%-8.5f,%-8.5f,%-8.5f,%-8.5f,%-8.5f,%-8.5f,%-8.5f\n'%\
				(x0,y0,B[2]*1000, B[5]*1000,emax*10**9,emin*10**9,thetastrain*180.0/np.pi,rotation*10**9,dilatacion*10**9,maxcizalla*10**9,maxgeo*10**9,np.abs(emax*10**9),np.abs(emin*10**9),90+thetastrain*180.0/np.pi,180+thetastrain*180.0/np.pi,270+thetastrain*180.0/np.pi))

    fidout.close()
    print 'fin'
