#!/usr/bin/env python

# coding:utf-8
from __future__ import print_function
import math
import numpy as np
#import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
import scipy as scipy
import scipy.integrate as integrate
import argparse
import time
import copy
#import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(description='Exact diagonalization of S=1/2 TFI in 1D')
    parser.add_argument('-N',metavar='N',dest='N',type=int,default=14,help='set Nsize (should be >=4)')
    parser.add_argument('-hi',metavar='hi',dest='hi',type=np.float,default=2.0,help='set hi=field_i/Hc (default hi=2.0, field_i=2.0*3.04438)')
    parser.add_argument('-hf',metavar='hf',dest='hf',type=np.float,default=0.0,help='set hf=field_f/Hc (default hf=0.0, field_f=0.0)')
#    parser.add_argument('-tau',metavar='tau',dest='tau',type=np.float,default=32.0,help='set tau (max time, default 32.0')
    return parser.parse_args()

def make_spin():
    S0 = scipy.sparse.csr_matrix(np.array([[1.0,0.0],[0.0,1.0]],dtype=np.float))
    Sx = scipy.sparse.csr_matrix(np.array([[0.0,1.0],[1.0,0.0]],dtype=np.float))
    Sy = scipy.sparse.csr_matrix(np.array([[0.0,-1.0j],[1.0j,0.0]],dtype=np.complex))
    Sz = scipy.sparse.csr_matrix(np.array([[1.0,0.0],[0.0,-1.0]],dtype=np.float))
    return S0,Sx,Sy,Sz

def make_list_2(N):
    list_site1 = [i for i in range(N)]
    list_site2 = [(i+1)%N for i in range(N)]
    list_Jzz = np.ones(N,dtype=np.float)
    return list_site1, list_site2, list_Jzz

def make_list_1(N):
    list_Jz = np.ones(N,dtype=np.float)
    return list_Jz

def make_list_1_stag(N):
    list_Jz = np.array([1.0 if i%2==0 else -1.0 for i in range(N)],dtype=np.float)
    return list_Jz

def make_ham_2(S0,Sz,N,Nbond,list_site1,list_site2,list_Jzz):
    Ham = scipy.sparse.csr_matrix((2**N,2**N),dtype=np.float)
    for bond in range(Nbond):
        i1 = list_site1[bond]
        i2 = list_site2[bond]
        Jzz = list_Jzz[bond]
        SzSz = 1
        for site in range(N):
            if site==i1 or site==i2:
                SzSz = scipy.sparse.kron(SzSz,Sz,format='csr')
            else:
                SzSz = scipy.sparse.kron(SzSz,S0,format='csr')
        Ham -= Jzz * SzSz
    return Ham

def make_ham_1(S0,Sz,N,list_Jz):
    Ham = scipy.sparse.csr_matrix((2**N,2**N),dtype=np.float)
    for site1 in range(N):
        ISz = 1
        Jz = list_Jz[site1]
        for site2 in range(N):
            if site2==site1:
                ISz = scipy.sparse.kron(ISz,Sz,format='csr')
            else:
                ISz = scipy.sparse.kron(ISz,S0,format='csr')
        Ham -= Jz * ISz
    return Ham

#----

def main():
    np.set_printoptions(threshold=10000)
    args = parse_args()
    N = args.N
    Nbond = N
    hi = args.hi
    hf = args.hf
#    tau = args.tau
#
    Hc = 1.0 ## 1D
#    Hc = 3.04438 ## 2D
#    tau = 32.0 # total time
#    dt = tau/1000.0
#
#    field_i = 2.0*Hc
#    field_f = 0.0
    field_i = hi*Hc
    field_f = hf*Hc
#    time_i = 0.0
#    time_f = tau
#
#    Nsteps = int(tau/dt+0.1)+1
#    Nsteps = 10+1
    Nsteps = 100+1
#    time_steps = [time_i+i*(time_f-time_i)/(Nsteps-1) for i in range(Nsteps)]
#    time_steps_for_U = [time_i+(i+0.5)*(time_f-time_i)/(Nsteps-1) for i in range(Nsteps)]
    field_steps = [field_i+i*(field_f-field_i)/(Nsteps-1) for i in range(Nsteps)]
#    field_steps_for_U = [field_i+(i+0.5)*(field_f-field_i)/(Nsteps-1) for i in range(Nsteps)]
#
    eps_break_Z2 = 0.0
#    eps_break_Z2 = 1e-4
#
    print("N=",N)
    print("Hc=",Hc)
    print("field_i=",field_i)
    print("field_f=",field_f)
#    print("tau=",tau)
    print("tau=infty")
#    print("dt=",dt)
    print("Nsteps=",Nsteps)
#    print("time_i=",time_i)
#    print("time_f=",time_f)
#    print("time_steps: t=",time_steps)
#    print("time_steps_for_U: t+dt/2=",time_steps_for_U)
    print("field_steps: h(t)=",field_steps)
#    print("field_steps_for_U: h(t+dt/2)=",field_steps_for_U)
    print("small -h_z\sigma^z: hz=",eps_break_Z2)

## prepare interaction
    start = time.time()
    S0, Sx, Sy, Sz = make_spin()
#    print(S0, Sx, Sy, Sz)
    list_site1, list_site2, list_Jzz = make_list_2(N)
    _, _, list_Jxx = make_list_2(N)
    list_Jx = make_list_1(N)
    list_Jz = make_list_1(N)
    list_Jx_stag = make_list_1_stag(N)
    list_Jz_stag = make_list_1_stag(N)
    print("site1",list_site1)
    print("site2",list_site2)
    print("Jzz",list_Jzz)
    print("Jx",list_Jx)
    end = time.time()
    print("time: prepare interaction",end - start)

## prepare Hamiltonian
    start = time.time()
    Ham_zz = make_ham_2(S0,Sz,N,Nbond,list_site1,list_site2,list_Jzz) ## FM Ising
#    Ham_zz = make_ham_2(S0,Sz,N,Nbond,list_site1,list_site2,-list_Jzz) ## AF Ising
    Ham_x = make_ham_1(S0,Sx,N,list_Jx)
    Op_xx = make_ham_2(S0,Sx,N,Nbond,list_site1,list_site2,list_Jxx)
    Op_zz = make_ham_2(S0,Sz,N,Nbond,list_site1,list_site2,list_Jzz)
    Op_x = make_ham_1(S0,Sx,N,list_Jx)
    Op_z = make_ham_1(S0,Sz,N,list_Jz)
    Op_x_stag = make_ham_1(S0,Sx,N,list_Jx_stag)
    Op_z_stag = make_ham_1(S0,Sz,N,list_Jz_stag)
    Op_Mz0Mz1 = copy.deepcopy(-Op_zz/Nbond) ## n.n. correlation <mz(0).mz(1)>
    Op_Mx0Mx1 = copy.deepcopy(-Op_xx/Nbond) ## n.n. correlation <mx(0).mx(1)>
    Op_Mx = copy.deepcopy(-Op_x/N)
    Op_Mz = copy.deepcopy(-Op_z/N)
    Op_Mx_stag = copy.deepcopy(-Op_x_stag/N)
    Op_Mz_stag = copy.deepcopy(-Op_z_stag/N)
    end = time.time()
    print("time: prepare Hamiltonian",end - start)

## prepare initial state: GS at t=0
    start = time.time()
    if np.abs(field_i)<1e-10:
        Ham0 = copy.deepcopy(Op_z) ## GS is |+z>
    else:
#        Ham0 = copy.deepcopy(Ham_zz + field_i*Ham_x) ## GS of TFI
        Ham0 = copy.deepcopy(Ham_zz + field_i*Ham_x + eps_break_Z2*Op_z) ## GS of TFI, break Z2 by hand
#
#    Ham0 = copy.deepcopy(Ham_zz + field_i*Ham_x) ## GS of TFI
#    Ham0 = copy.deepcopy(Op_x) ## GS is |+x>
#    Ham0 = copy.deepcopy(Op_z) ## GS is |+z>
#    Ham0 = copy.deepcopy(Op_x_stag) ## GS is staggered AF along x axis
#    Ham0 = copy.deepcopy(Op_z_stag) ## GS is staggered AF along z axis
#
    ene,vec = scipy.sparse.linalg.eigsh(Ham0,which='SA',k=2) # real Ham
    print("energy(t=0)",ene[0],ene[1])
#    print("vector:",vec[:,0],vec[:,1])
    end = time.time()
    print("time: prepare intial state",end - start)

## calculate dynamics
    start = time.time()
    psi = copy.deepcopy(vec[:,0])
    ret = []
    ret.append(psi)
#    for i in range(Nsteps):
    for i in range(1,Nsteps):
#        time_m = time_steps[i]
        field_m = field_steps[i]
#        Ham = copy.deepcopy(Ham_zz + field_m*Ham_x)
        Ham = copy.deepcopy(Ham_zz + field_m*Ham_x + eps_break_Z2*Op_z) ## break Z2 by hand
        ene,vec = scipy.sparse.linalg.eigsh(Ham,which='SA',k=2)
        psi = vec[:,0]
        ret.append(psi)
#        print("psi",i,psi)
#    print(ret)
    ene0 = [(np.dot(np.conjugate(ret[i]),Ham0.dot(ret[i])) / np.dot(np.conjugate(ret[i]),ret[i])).real for i in range(Nsteps)]
    ene1 = [(np.dot(np.conjugate(ret[i]),(Ham_zz+field_steps[Nsteps-1]*Ham_x).dot(ret[i])) / np.dot(np.conjugate(ret[i]),ret[i])).real for i in range(Nsteps)]
#    enezz = [(np.dot(np.conjugate(ret[i]),Ham_zz.dot(ret[i])) / np.dot(np.conjugate(ret[i]),ret[i])).real for i in range(Nsteps)]
#    enex = [(np.dot(np.conjugate(ret[i]),Ham_x.dot(ret[i])) / np.dot(np.conjugate(ret[i]),ret[i])).real for i in range(Nsteps)]
    ene = [(np.dot(np.conjugate(ret[i]),(Ham_zz+field_steps[i]*Ham_x).dot(ret[i])) / np.dot(np.conjugate(ret[i]),ret[i])).real for i in range(Nsteps)]
    mx0mx1 = [(np.dot(np.conjugate(ret[i]),Op_Mx0Mx1.dot(ret[i])) / np.dot(np.conjugate(ret[i]),ret[i])).real for i in range(Nsteps)]
    mz0mz1 = [(np.dot(np.conjugate(ret[i]),Op_Mz0Mz1.dot(ret[i])) / np.dot(np.conjugate(ret[i]),ret[i])).real for i in range(Nsteps)]
    mx = [(np.dot(np.conjugate(ret[i]),Op_Mx.dot(ret[i])) / np.dot(np.conjugate(ret[i]),ret[i])).real for i in range(Nsteps)]
    mz = [(np.dot(np.conjugate(ret[i]),Op_Mz.dot(ret[i])) / np.dot(np.conjugate(ret[i]),ret[i])).real for i in range(Nsteps)]
    mx_stag = [(np.dot(np.conjugate(ret[i]),Op_Mx_stag.dot(ret[i])) / np.dot(np.conjugate(ret[i]),ret[i])).real for i in range(Nsteps)]
    mz_stag = [(np.dot(np.conjugate(ret[i]),Op_Mz_stag.dot(ret[i])) / np.dot(np.conjugate(ret[i]),ret[i])).real for i in range(Nsteps)]
    loschmidt_echo = [(np.abs(np.dot(np.conjugate(ret[i]),ret[0]))**2 / np.dot(np.conjugate(ret[i]),ret[i])*np.dot(np.conjugate(ret[0]),ret[0])).real for i in range(Nsteps)]
    print("ene0",ene0)
    print("ene",ene)
    print("mz0mz1",mz0mz1)
    print("mx0mx1",mx0mx1)
    print("mx",mx)
    print("mz",mz)
    print("mx_stag",mx_stag)
    print("mz_stag",mz_stag)
    print("loschmidt_echo",loschmidt_echo)
#    print("kink density of the final state: N tau KinkDens",N,tau,0.5*(1.0-mz0mz1[Nsteps-1]))
    print("kink density of the final state: N tau KinkDens",N,"inf",0.5*(1.0-mz0mz1[Nsteps-1]))
    end = time.time()
    print("time: calculate dynamics",end - start)

## plot energy evolution
    start = time.time()
#
    fig0 = plt.figure()
    fig0.suptitle("energy0")
    plt.plot(field_steps,ene0)
    plt.xlabel("field")
    plt.gca().invert_xaxis()
    fig0.savefig("fig_ene0.png")
#
    fig11 = plt.figure()
    fig11.suptitle("energy1")
    plt.plot(field_steps,ene1)
    plt.xlabel("field")
    plt.gca().invert_xaxis()
    fig11.savefig("fig_ene1.png")
#
    fig1 = plt.figure()
    fig1.suptitle("energy")
    plt.plot(field_steps,ene)
    plt.xlabel("field")
    plt.gca().invert_xaxis()
    fig1.savefig("fig_ene.png")
#
    fig2 = plt.figure()
    fig2.suptitle("mx0mx1")
    plt.plot(field_steps,mx0mx1)
    plt.xlabel("field")
    plt.gca().invert_xaxis()
    fig2.savefig("fig_mx0mx1.png")
#
    fig3 = plt.figure()
    fig3.suptitle("mz0mz1")
    plt.plot(field_steps,mz0mz1)
    plt.xlabel("field")
    plt.gca().invert_xaxis()
    fig3.savefig("fig_mz0mz1.png")
#
    fig4 = plt.figure()
    fig4.suptitle("mx")
    plt.plot(field_steps,mx)
    plt.xlabel("field")
    plt.gca().invert_xaxis()
    fig4.savefig("fig_mx.png")
#
    fig5 = plt.figure()
    fig5.suptitle("mz")
#    fig5.suptitle("mz (small $-h_z\sigma^z$ is added)")
    plt.plot(field_steps,mz)
    plt.xlabel("field")
    plt.gca().invert_xaxis()
    fig5.savefig("fig_mz.png")
#
#    fig6 = plt.figure()
#    fig6.suptitle("mx_stag")
#    plt.plot(field_steps,mx_stag)
#    plt.xlabel("field")
#    plt.gca().invert_xaxis()
#    fig6.savefig("fig_mx_stag.png")
#
#    fig7 = plt.figure()
#    fig7.suptitle("mz_stag")
#    plt.plot(field_steps,mz_stag)
#    plt.xlabel("field")
#    plt.gca().invert_xaxis()
#    fig7.savefig("fig_mz_stag.png")
#
    fig8 = plt.figure()
    fig8.suptitle("Loschmidt echo")
    plt.plot(field_steps,loschmidt_echo)
    plt.xlabel("field")
    plt.gca().invert_xaxis()
    fig8.savefig("fig_loschmidt_echo.png")
#
    fig9 = plt.figure()
    fig9.suptitle("logarithmic Loschmidt echo")
    plt.plot(field_steps,-np.log(loschmidt_echo)/N)
    plt.xlabel("field")
    plt.gca().invert_xaxis()
    fig9.savefig("fig_log_loschmidt_echo.png")
#
    fig10 = plt.figure()
    fig10.suptitle("kink density")
    plt.plot(field_steps,0.5*(np.ones(Nsteps,dtype=np.float)-mz0mz1))
    plt.xlabel("field")
    plt.gca().invert_xaxis()
    fig10.savefig("fig_kink_density.png")
#
#    plt.show()
    end = time.time()
    print("time: plot",end - start)

if __name__ == "__main__":
    main()
