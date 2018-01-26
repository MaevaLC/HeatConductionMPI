# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 16:20:49 2018

@author: m.leclech
"""

import matplotlib.pyplot as plt

fileFTCS1 = open("FTCS1-0.5.txt", 'r')
fileFTCS2 = open("FTCS1-0.05.txt", 'r')
fileFTCS3 = open("FTCS1-0.005.txt", 'r')

fileAN5 = open("Analytical-0.500000.txt", 'r')

dataFTCS1 = fileFTCS1.readlines()
dataFTCS2 = fileFTCS2.readlines()
dataFTCS3 = fileFTCS3.readlines()

dataAN5 = fileAN5.readlines()

x1 = []
x2 = []
x3 = []

FTCS1 = []
FTCS2 = []
FTCS3 = []

AN5 = []

i = 0
while ( i < len(dataFTCS1) ):
    data_FTCS1_splited = dataFTCS1[i].split(" ")    
    x1.append(float(data_FTCS1_splited[0]))    
    FTCS1.append( float(data_FTCS1_splited[1]) )
    i+=1

j = 0
while ( j < len(dataFTCS2) ):
    data_FTCS2_splited = dataFTCS2[j].split(" ")
    data_AN5_splited = dataAN5[j].split(" ")
    x2.append(float(data_FTCS2_splited[0]))    
    FTCS2.append( float(data_FTCS2_splited[1]) )
    AN5.append( float(data_AN5_splited[1]) )
    j+=1
    
k = 0
while ( k < len(dataFTCS3) ):
    data_FTCS3_splited = dataFTCS3[k].split(" ")    
    x3.append(float(data_FTCS3_splited[0]))    
    FTCS3.append( float(data_FTCS3_splited[1]) )
    k+=1

plt.subplot(111)    
plt.plot(x1, FTCS1, 'r', label='dx = 0.5')
#plt.plot(x2, FTCS2, 'g', label='dx = 0.05')
#plt.plot(x3, FTCS3, 'b', label='dx = 0.005')
plt.plot(x2, AN5, 'm', label='analytical')
plt.title("FTCS")
plt.xlim((0, 1))
#plt.ylim((100, 300))
plt.legend(loc ='upper center')
plt.show()

fileFTCS1.close()
fileFTCS2.close()
fileFTCS3.close()
fileAN5.close()