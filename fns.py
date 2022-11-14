import numpy as np
from numpy import copy
import para
from compressible import Compressible
from derivative import *
from boundary import *
import h5py


def update_primitive(compress):

        compress.rho = compress.Q[0]

        compress.ux = compress.Q[1]/compress.rho
        compress.uz = compress.Q[2]/compress.rho
        imposeBC_u(compress)

        compress.e_t = compress.Q[3]/compress.rho

        compress.T = (compress.e_t - ((1/2)*(compress.ux**2 + \
                    compress.uz**2)) - para.g*para.Z_mesh)/para.Cv
        imposeBC_T(compress)

        compress.pressure()
        imposeBC_p(compress)
        
        compress.density()
        compress.specific_energy()

        pass

def time_advance_single_step(dt, compress=Compressible()):

    compress.compute_rhs()
    compress.Q += compress.rhs*dt

    update_primitive(compress)

    return compress

def time_advance_euler(compress=Compressible()):

    t = para.tinit
    j = 0

    E1 = np.array([])
    E2 = np.array([])
    E3 = np.array([])

    N = np.array([])

    while t <= para.tfinal:

        # print(t, total_energy(compress))

        # Plots
        if ((para.tstep[j]-t)/para.dt) <= para.dt:

            hf = h5py.File("output/2D_%d.h5" % (j), 'w')
            hf.create_dataset('p', data=compress.p)
            hf.create_dataset('ux', data=compress.ux)
            hf.create_dataset('uz', data=compress.uz)
            hf.create_dataset('T', data=compress.T)
            hf.close()

            j=j+1

            pass

        compress=time_advance_single_step(para.dt, compress)

        E1 = np.append(E1, total_energy(compress))
        E2 = np.append(E2, total_kinetic_energy(compress))
        E3 = np.append(E3, total_internal_energy(compress))

        N = np.append(N, Nusselt_number(compress))

        t = np.round(t+para.dt, para.n1)

        if np.sum(E1) > 1e16:

            print ('Try different parameters..Energy became infinity, code blew up at t = ', t)    

            break

        pass

    hf2 = h5py.File("output\E&Nu.h5", 'w')
    hf2.create_dataset('E1', data=E1)
    hf2.create_dataset('E2', data=E2)
    hf2.create_dataset('E3', data=E3)
    hf2.create_dataset('N', data=N)

    pass


def time_advance_RK2(compress=Compressible()):

    t = para.tinit
    j = 0

    E1 = np.array([])
    E2 = np.array([])
    E3 = np.array([])

    N = np.array([])

    while t <= para.tfinal:

        # print(t, total_energy(compress))

        # Plots
        if ((para.tstep[j]-t)/para.dt) <= para.dt:

            hf = h5py.File("output/2D_%d.h5" % (j), 'w')
            hf.create_dataset('p', data=compress.p)
            hf.create_dataset('ux', data=compress.ux)
            hf.create_dataset('uz', data=compress.uz)
            hf.create_dataset('T', data=compress.T)
            hf.close()

            j=j+1

            pass

        compress.compute_rhs()
        compress.Q_copy = copy(compress.Q)
        compress.Q += compress.rhs*para.dt/2
        update_primitive(compress)

        compress.compute_rhs()       
        compress.Q = compress.Q_copy + compress.rhs*para.dt    
        update_primitive(compress)

        E1 = np.append(E1, total_energy(compress))
        E2 = np.append(E2, total_kinetic_energy(compress))
        E3 = np.append(E3, total_internal_energy(compress))

        N = np.append(N, Nusselt_number(compress))

        t = np.round(t+para.dt, para.n1)

        if np.sum(E1) > 1e16:

            print ('Try different parameters..Energy became infinity, code blew up at t = ', t)    

            break

        pass

    hf2 = h5py.File("output/E&Nu.h5", 'w')
    hf2.create_dataset('E1', data=E1)
    hf2.create_dataset('E2', data=E2)
    hf2.create_dataset('E3', data=E3)
    hf2.create_dataset('N', data=N)

    pass


def total_energy(compress=Compressible()):
    # Total energy integral

    return np.sum(compress.rho*compress.e_t)/(para.Nx*para.Nz)


def total_kinetic_energy(compress=Compressible()):
    # Total kinetic energy integral

    return np.sum(compress.rho*(1/2)*(para.A**2*compress.ux**2 + compress.uz**2))/(para.Nx*para.Nz)

def total_internal_energy(compress=Compressible()):
    # Total kinetic energy integral

    return np.sum(compress.rho*(para.Cv*compress.T))/(para.Nx*para.Nz)


def Nusselt_number(compress=Compressible()):

    return - np.sum(dfz_c(compress.T)[:,0])/(para.Nx)