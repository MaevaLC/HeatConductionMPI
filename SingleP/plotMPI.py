# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 16:20:49 2018

@author: m.leclech
"""

import matplotlib.pyplot as plt

fileFTCS1 = open("FTCS-0.5.txt", 'r')
fileFTCS2 = open("FTCS-0.05.txt", 'r')
fileFTCS3 = open("FTCS-0.005.txt", 'r')

fileLAAS1 = open("LAAS-0.5.txt", 'r')
fileLAAS2 = open("LAAS-0.05.txt", 'r')
fileLAAS3 = open("LAAS-0.005.txt", 'r')

fileCRANK1 = open("CRANK-0.5.txt", 'r')
fileCRANK2 = open("CRANK-0.05.txt", 'r')
fileCRANK3 = open("CRANK-0.005.txt", 'r')

fileAN5 = open("Analytical-0.500000.txt", 'r')

dataFTCS1 = fileFTCS1.readlines()
dataFTCS2 = fileFTCS2.readlines()
dataFTCS3 = fileFTCS3.readlines()

dataLAAS1 = fileLAAS1.readlines()
dataLAAS2 = fileLAAS2.readlines()
dataLAAS3 = fileLAAS3.readlines()

dataCRANK1 = fileCRANK1.readlines()
dataCRANK2 = fileCRANK2.readlines()
dataCRANK3 = fileCRANK3.readlines()

dataAN5 = fileAN5.readlines()

x1 = []
x2 = []
x3 = []

FTCS1 = []
FTCS2 = []
FTCS3 = []

LAAS1 = []
LAAS2 = []
LAAS3 = []

CRANK1 = []
CRANK2 = []
CRANK3 = []

AN5 = []

i = 0
while ( i < len(dataFTCS1) ):
    data_FTCS1_splited = dataFTCS1[i].split(" ")
    data_LAAS1_splited = dataLAAS1[i].split(" ")
    data_CRANK1_splited = dataCRANK1[i].split(" ")
    x1.append(float(data_FTCS1_splited[0]))    
    FTCS1.append( float(data_FTCS1_splited[1]) )
    LAAS1.append( float(data_LAAS1_splited[1]) )
    CRANK1.append( float(data_CRANK1_splited[1]) )
    i+=1

j = 0
while ( j < len(dataFTCS2) ):
    data_FTCS2_splited = dataFTCS2[j].split(" ")
    data_LAAS2_splited = dataLAAS2[j].split(" ")
    data_CRANK2_splited = dataCRANK2[j].split(" ")
    data_AN5_splited = dataAN5[j].split(" ")
    x2.append(float(data_FTCS2_splited[0]))    
    FTCS2.append( float(data_FTCS2_splited[1]) )
    LAAS2.append( float(data_LAAS2_splited[1]) )
    CRANK2.append( float(data_CRANK2_splited[1]) )
    AN5.append( float(data_AN5_splited[1]) )
    j+=1
    
k = 0
while ( k < len(dataFTCS3) ):
    data_FTCS3_splited = dataFTCS3[k].split(" ")
    data_LAAS3_splited = dataLAAS3[k].split(" ")
    data_CRANK3_splited = dataCRANK3[k].split(" ")
    x3.append(float(data_FTCS3_splited[0]))    
    FTCS3.append( float(data_FTCS3_splited[1]) )
    LAAS3.append( float(data_LAAS3_splited[1]) )
    CRANK3.append( float(data_CRANK3_splited[1]) )
    k+=1

plt.subplot(311)    
plt.plot(x1, FTCS1, 'r', label='dx = 0.5')
plt.plot(x3, FTCS3, 'g', label='dx = 0.005')
plt.plot(x2, AN5, 'm', label='analytical')
plt.title("FTCS")
plt.xlim((0, 1))
#plt.ylim((100, 300))
plt.legend(loc ='upper center')


plt.subplot(312)    
plt.plot(x1, LAAS1, 'r', label='dx = 0.5')
plt.plot(x3, LAAS3, 'g', label='dx = 0.005')
plt.plot(x2, AN5, 'm', label='analytical')
plt.title("LAAS")
plt.xlim((0, 1))
plt.ylim((140, 300))
plt.legend(loc ='upper center')


plt.subplot(313)    
plt.plot(x1, CRANK1, 'r', label='dx = 0.5')
plt.plot(x3, CRANK3, 'g', label='dx = 0.005')
plt.plot(x2, AN5, 'm', label='analytical')
plt.title("CRANK")
plt.xlim((0, 1))
plt.ylim((140, 300))
plt.legend(loc ='upper center')
plt.show()

plt.subplot(211)    
plt.plot(x1, FTCS1, 'r', label='FTCS')
plt.plot(x1, LAAS1, 'g', label='LAAS')
plt.plot(x1, CRANK1, 'b', label='CRANK')
plt.plot(x2, AN5, 'm', label='analytical')
plt.title("dx = 0.5")
plt.xlim((0, 1))
plt.ylim((140, 300))
plt.legend(loc ='upper center')


plt.subplot(212)    
plt.plot(x3, FTCS3, 'r', label='FTCS')
plt.plot(x3, LAAS3, 'g', label='LAAS')
plt.plot(x3, CRANK3, 'b', label='CRANK')
plt.plot(x2, AN5, 'm', label='analytical')
plt.title("dx = 0.005")
plt.xlim((0, 1))
plt.ylim((140, 300))
plt.legend(loc ='upper center')
plt.show()

fileFTCS1.close()
fileFTCS2.close()
fileFTCS3.close()
fileLAAS1.close()
fileLAAS2.close()
fileLAAS3.close()
fileCRANK1.close()
fileCRANK2.close()
fileCRANK3.close()
fileAN5.close()