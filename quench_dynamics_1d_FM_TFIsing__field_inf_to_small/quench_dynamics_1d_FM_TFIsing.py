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
    parser = argparse.ArgumentParser(description='Dynamics of S=1/2 TFI chain')
    parser.add_argument('-N',metavar='N',dest='N',type=int,default=10,help='set Nsize (should be >=4)')
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
# analytical result
#
def epsion(J,h,k):
    return 2.0*J*np.sqrt(1.0+h**2-2.0*h*np.cos(k))

def diff_m_field(J,h0,h1,k,t):
    e0 = epsion(J,h0,k)
    e1 = epsion(J,h1,k)
    denom = e1**2 * e0
    numer0 = (h1-h0) * (np.sin(k))**2 * np.cos(2.0*e1*t)
    numer1 = (h1*h0+1.0-(h1+h0)*np.cos(k)) * (h1-np.cos(k))
    return (numer1-numer0)/denom * (4.0/np.pi)

def diff_m_field_main(k,c):
    return diff_m_field(c[0],c[1],c[2],k,c[3])

def analytical_m_field(J,h0,h1,tau,dt):
    list_vt = []
    list_m_field = []
    vmax = 2.0*J*np.min([1.0,h1])
    tmax = tau/vmax
    N = int(tmax/dt+0.1)+1
    for steps in range(N):
        t = steps*dt
        m_field = integrate.quad(lambda k,c : \
            diff_m_field_main(k,c), 0.0, 2.0*np.pi, args=[J,h0,h1,t])
#        print("{0:.16f} {1:.16f} {2:.16f} {3:.16f}".format(t,vmax*t,m_field[0],m_field[1]))
        list_vt.append(vmax*t)
        list_m_field.append(m_field[0])
    return list_vt, list_m_field

def diff_m_ising(J,h0,h1,k):
    e0 = epsion(J,h0,k)
    e1 = epsion(J,h1,k)
    e1prime = 4.0*J**2*h1*np.sin(k)/e1
    cosDk = 4.0*J**2*(h1*h0-(h1+h0)*np.cos(k)+1.0)/e1/e0
    return e1prime*np.log(np.abs(cosDk))/np.pi

def diff_m_ising_main(k,c):
    return diff_m_ising(c[0],c[1],c[2],k)

def calc_CFF(h0,h1):
    assert np.abs(h0)<=1.0 and np.abs(h1)<=1.0, "!!! choose |h0|,|h1|<1"
    numer = 1.0-h1*h0+np.sqrt((1.0-h1**2)*(1.0-h0**2))
    denom = 2.0*np.sqrt((1.0-h1*h0)*np.sqrt(1.0-h0**2))
    return np.sqrt(numer/denom)

def analytical_m_ising(J,h0,h1,tau,dt):
    list_vt = []
    list_m_ising = []
    vmax = 2.0*J*np.min([1.0,h1])
    tmax = tau/vmax
    N = int(tmax/dt+0.1)+1
    CFF = calc_CFF(h0,h1)
    integral_m_ising = integrate.quad(lambda k,c : \
        diff_m_ising_main(k,c), 0.0, np.pi, args=[J,h0,h1])
#    print("# m(t)~a*exp(bt)")
#    print("# a~ {0:.16f}".format(CFF))
#    print("# b~ {0:.16f} {1:.16f}".format(integral_m_ising[0],integral_m_ising[1]))
    for steps in range(N):
        t = steps*dt
        m_ising = CFF*np.exp(t*integral_m_ising[0])
#        print("{0:.16f} {1:.16f} {2:.16f}".format(t,vmax*t,m_ising))
        list_vt.append(vmax*t)
        list_m_ising.append(m_ising)
    return list_vt, list_m_ising

def diff_log_loschmidt_echo(J,h0,h1,k,t):
    e0 = epsion(J,h0,k)
    e1 = epsion(J,h1,k)
    denom = e1 * e0
    numer = 4.0*J**2 * (h1-h0) * np.sin(k) * np.sin(e1*t)
    return np.log(1.0 - (numer/denom)**2) / (-2.0*np.pi)

def diff_log_loschmidt_echo_main(k,c):
    return diff_log_loschmidt_echo(c[0],c[1],c[2],k,c[3])

def analytical_log_loschmidt_echo(J,h0,h1,tau,dt):
    list_vt = []
    list_log_loschmidt_echo = []
    vmax = 2.0*J*np.min([1.0,h1])
    tmax = tau/vmax
    N = int(tmax/dt+0.1)+1
    for steps in range(N):
        t = steps*dt
        log_loschmidt_echo = integrate.quad(lambda k,c : \
            diff_log_loschmidt_echo_main(k,c), 0.0, np.pi, args=[J,h0,h1,t])
#        print("{0:.16f} {1:.16f} {2:.16f} {3:.16f}".format(t,vmax*t,log_loschmidt_echo[0],log_loschmidt_echo[1]))
        list_vt.append(vmax*t)
        list_log_loschmidt_echo.append(log_loschmidt_echo[0])
    return list_vt, list_log_loschmidt_echo
#----

def main():
    np.set_printoptions(threshold=10000)
    args = parse_args()
    N = args.N
    Nbond = N
#
    field_i = 1e10
    field_f = 0.1
    print("N=",N)
    print("field_i=",field_i)
    print("field_f=",field_f)
#
    J = 1.0
    tau = 1.0
    dt = 0.01
    print("J=",J)
    print("tau=",tau)
    print("dt=",dt)

## analytical
    start = time.time()
    list_vt, list_m_field = analytical_m_field(J,field_i,field_f,tau,dt)
#    if np.abs(field_i)<=1.0 and np.abs(field_f)<=1.0:
    if np.abs(field_i)<=1.0 and np.abs(field_f)<=1.0 and np.abs(field_i)<1e-10: ## when an initial state explicitly breaks Z2 symmetry
        _, list_m_ising = analytical_m_ising(J,field_i,field_f,tau,dt)
        _, list_log_loschmidt_echo = analytical_log_loschmidt_echo(J,field_i,field_f,tau,dt)
    end = time.time()
    print("time: analytical",end - start)

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
        Ham0 = copy.deepcopy(Ham_zz + field_i*Ham_x) ## GS of TFI
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
#
    h1 = field_f
    vmax = 2.0*J*np.min([1.0,h1])
    tmax = tau/vmax
    Nsteps = int(tmax/dt+0.1)+1
#
    timei = 0.0
    timef = tmax
    steps = Nsteps
    Ham = copy.deepcopy(Ham_zz + field_f*Ham_x)
#    print(Ham)
    psi0 = copy.deepcopy(vec[:,0])
    ret = scipy.sparse.linalg.expm_multiply((-1j)*Ham,vec[:,0],start=timei,stop=timef,num=steps,endpoint=True)
#    print(ret)
    ene0 = [(np.dot(np.conjugate(ret[i]),Ham0.dot(ret[i])) / np.dot(np.conjugate(ret[i]),ret[i])).real for i in range(steps)]
    ene = [(np.dot(np.conjugate(ret[i]),Ham.dot(ret[i])) / np.dot(np.conjugate(ret[i]),ret[i])).real for i in range(steps)]
    mx0mx1 = [(np.dot(np.conjugate(ret[i]),Op_Mx0Mx1.dot(ret[i])) / np.dot(np.conjugate(ret[i]),ret[i])).real for i in range(steps)]
    mz0mz1 = [(np.dot(np.conjugate(ret[i]),Op_Mz0Mz1.dot(ret[i])) / np.dot(np.conjugate(ret[i]),ret[i])).real for i in range(steps)]
    mx = [(np.dot(np.conjugate(ret[i]),Op_Mx.dot(ret[i])) / np.dot(np.conjugate(ret[i]),ret[i])).real for i in range(steps)]
    mz = [(np.dot(np.conjugate(ret[i]),Op_Mz.dot(ret[i])) / np.dot(np.conjugate(ret[i]),ret[i])).real for i in range(steps)]
    mx_stag = [(np.dot(np.conjugate(ret[i]),Op_Mx_stag.dot(ret[i])) / np.dot(np.conjugate(ret[i]),ret[i])).real for i in range(steps)]
    mz_stag = [(np.dot(np.conjugate(ret[i]),Op_Mz_stag.dot(ret[i])) / np.dot(np.conjugate(ret[i]),ret[i])).real for i in range(steps)]
    loschmidt_echo = [(np.abs(np.dot(np.conjugate(ret[i]),psi0))**2 / np.dot(np.conjugate(ret[i]),ret[i])*np.dot(np.conjugate(psi0),psi0)).real for i in range(steps)]
    print("ene0",ene0)
    print("ene",ene)
    print("mz0mz1",mz0mz1)
    print("mx0mx1",mx0mx1)
    print("mx",mx)
    print("mz",mz)
    print("mx_stag",mx_stag)
    print("mz_stag",mz_stag)
    print("loschmidt_echo",loschmidt_echo)
    end = time.time()
    print("time: calculate dynamics",end - start)

## plot energy evolution
    start = time.time()
#    time_steps = [timei+i*(timef-timei)/(steps-1) for i in range(steps)]
    time_steps = [vmax*(timei+i*(timef-timei)/(steps-1)) for i in range(steps)]
#
#    fig0 = plt.figure()
#    fig0.suptitle("energy0")
#    plt.plot(time_steps,ene0)
#    plt.xlabel("$v_{max}t$")
#    fig0.savefig("fig_ene0.png")
#
#    fig1 = plt.figure()
#    fig1.suptitle("energy")
#    plt.plot(time_steps,ene)
#    plt.xlabel("$v_{max}t$")
#    fig1.savefig("fig_ene.png")
#
#    fig2 = plt.figure()
#    fig2.suptitle("mx0mx1")
#    plt.plot(time_steps,mx0mx1)
#    plt.xlabel("$v_{max}t$")
#    fig2.savefig("fig_mx0mx1.png")
#
    fig102 = plt.figure()
    fig102.suptitle("mx0mx1")
    plt.plot(time_steps/vmax/np.pi*4.0,mx0mx1)
    plt.xlabel("$v_{max}t$")
    fig102.savefig("fig_mx0mx1_vs_t.png")
#
#    fig3 = plt.figure()
#    fig3.suptitle("mz0mz1")
#    plt.plot(time_steps,mz0mz1)
#    plt.xlabel("$v_{max}t$")
#    fig3.savefig("fig_mz0mz1.png")
#
    fig103 = plt.figure()
    fig103.suptitle("mz0mz1/$h_f$")
    plt.plot(time_steps/vmax/np.pi*4.0,np.array(mz0mz1,dtype=np.float)/field_f)
    plt.xlabel("$t/\pi$")
    fig103.savefig("fig_mz0mz1_vs_t.png")
#
#    fig4 = plt.figure()
#    fig4.suptitle("mx")
#    plt.plot(time_steps,mx,label="ED")
#    plt.plot(time_steps,list_m_field,label="Analytical")
#    plt.legend(bbox_to_anchor=(1,1),loc='upper right',borderaxespad=1,fontsize=12)
#    plt.xlabel("$v_{max}t$")
#    fig4.savefig("fig_mx.png")
#
    fig104 = plt.figure()
    fig104.suptitle("mx")
    plt.plot(time_steps/vmax/np.pi*4.0,mx,label="ED")
    plt.plot(time_steps/vmax/np.pi*4.0,list_m_field,label="Analytical")
    plt.legend(bbox_to_anchor=(1,1),loc='upper right',borderaxespad=1,fontsize=12)
    plt.xlabel("$t/\pi$")
    fig104.savefig("fig_mx_vs_t.png")
#
#    fig5 = plt.figure()
#    fig5.suptitle("mz")
#    plt.plot(time_steps,mz,label="ED")
##    if np.abs(field_i)<=1.0 and np.abs(field_f)<=1.0:
#    if np.abs(field_i)<=1.0 and np.abs(field_f)<=1.0 and np.abs(field_i)<1e-10: ## when an initial state explicitly breaks Z2 symmetry
#        plt.plot(time_steps,list_m_ising,label="Asymptote")
#    plt.legend(bbox_to_anchor=(1,1),loc='upper right',borderaxespad=1,fontsize=12)
#    plt.xlabel("$v_{max}t$")
#    fig5.savefig("fig_mz.png")
#
#    fig6 = plt.figure()
#    fig6.suptitle("mx_stag")
#    plt.plot(time_steps,mx_stag)
#    plt.xlabel("$v_{max}t$")
#    fig6.savefig("fig_mx_stag.png")
#
#    fig7 = plt.figure()
#    fig7.suptitle("mz_stag")
#    plt.plot(time_steps,mz_stag)
#    plt.xlabel("$v_{max}t$")
#    fig7.savefig("fig_mz_stag.png")
#
#    fig8 = plt.figure()
#    fig8.suptitle("Loschmidt echo")
#    plt.plot(time_steps,loschmidt_echo)
#    plt.xlabel("$v_{max}t$")
#    fig8.savefig("fig_loschmidt_echo.png")
#
#    fig9 = plt.figure()
#    fig9.suptitle("logarithmic Loschmidt echo")
#    plt.plot(time_steps,-np.log(loschmidt_echo)/N,label="ED")
#    if np.abs(field_i)<1e-10: ## when an initial state explicitly breaks Z2 symmetry
#        plt.plot(time_steps,list_log_loschmidt_echo,label="Analytical")
#    plt.legend(bbox_to_anchor=(1,1),loc='upper right',borderaxespad=1,fontsize=12)
#    plt.xlabel("$v_{max}t$")
#    fig9.savefig("fig_loschmidt_echo.png")
#
    fig109 = plt.figure()
    fig109.suptitle("logarithmic Loschmidt echo")
    plt.plot(time_steps/vmax/np.pi*4.0,-np.log(loschmidt_echo)/N,label="ED")
    if np.abs(field_i)<1e-10: ## when an initial state explicitly breaks Z2 symmetry
        plt.plot(time_steps/vmax/np.pi*4.0,list_log_loschmidt_echo,label="Analytical")
    plt.legend(bbox_to_anchor=(1,1),loc='upper right',borderaxespad=1,fontsize=12)
    plt.xlabel("$t/\pi$")
    fig109.savefig("fig_loschmidt_echo_vs_t.png")
#
#    plt.show()
    end = time.time()
    print("time: plot",end - start)

if __name__ == "__main__":
    main()
