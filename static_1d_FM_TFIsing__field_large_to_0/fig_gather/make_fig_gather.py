#!/usr/bin/env python

# coding:utf-8
from __future__ import print_function
#import sys
import re
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

all_files = glob.glob('../*/dat_N*_tau_inf')

list_N = []
list_mx = []
list_mz0mz1 = []
list_ene = []
for file_name in all_files:
#    N = file_name.replace("dat_N","")
    N = re.sub(".*dat_N","",file_name)
    N = N.replace("_tau_inf","")
    list_N.append(int(N))
    print(file_name,N)
#    file = open(sys.argv[1])
#    file = open('dat_N10_tau_inf')
    file = open(file_name)
    lines = file.readlines()
    file.close()
    for line in lines:
        if line.startswith("mx ["):
            line_mx = line[:-1]
            line_mx = line_mx.replace("mx [","")
            line_mx = line_mx.replace("]","")
#            list_mx = np.fromstring(line_mx,dtype=np.float,sep=',')
            list_mx.append(np.fromstring(line_mx,dtype=np.float,sep=','))
        if line.startswith("mz0mz1 ["):
            line_mz0mz1 = line[:-1]
            line_mz0mz1 = line_mz0mz1.replace("mz0mz1 [","")
            line_mz0mz1 = line_mz0mz1.replace("]","")
            list_mz0mz1.append(np.fromstring(line_mz0mz1,dtype=np.float,sep=','))
        if line.startswith("ene ["):
            line_ene = line[:-1]
            line_ene = line_ene.replace("ene [","")
            line_ene = line_ene.replace("]","")
            list_ene.append(np.fromstring(line_ene,dtype=np.float,sep=','))
        if line.startswith("field_steps: h(t)= ["):
            line_h = line[:-1]
            line_h = line_h.replace("field_steps: h(t)= [","")
            line_h = line_h.replace("]","")
            list_h = np.fromstring(line_h,dtype=np.float,sep=',')

list_enedens = []
for i in range(len(list_N)):
    list_enedens.append([x/list_N[i] for x in list_ene[i]])

print("h",list_h)
for i in range(len(list_N)):
    print("N mx",list_N[i],list_mx[i])
    print("N mz0mz1",list_N[i],list_mz0mz1[i])
    print("N enedens",list_N[i],list_enedens[i])

fig0 = plt.figure()
fig0.suptitle("mx")
for i in range(len(list_N)):
    plt.plot(list_h,list_mx[i],label=list_N[i])
plt.xlabel("field")
plt.legend(bbox_to_anchor=(1,1),loc='upper right',borderaxespad=1)
plt.gca().invert_xaxis()
fig0.savefig("fig_mx.png")

fig1 = plt.figure()
fig1.suptitle("mz0mz1")
for i in range(len(list_N)):
    plt.plot(list_h,list_mz0mz1[i],label=list_N[i])
plt.xlabel("field")
plt.legend(bbox_to_anchor=(1,0),loc='lower right',borderaxespad=1)
plt.gca().invert_xaxis()
fig1.savefig("fig_mz0mz1.png")

fig2 = plt.figure()
fig2.suptitle("enedens")
for i in range(len(list_N)):
    plt.plot(list_h,list_enedens[i],label=list_N[i])
plt.xlabel("field")
plt.legend(bbox_to_anchor=(1,0),loc='lower right',borderaxespad=1)
plt.gca().invert_xaxis()
fig2.savefig("fig_enedens.png")
