import numpy as np
import para
from derivative import *
from numpy import copy
import h5py

class Compressible:

    def __init__(self):
        
        # Velocity vector components
        self.ux = []
        self.uz = []

        self.DuxDx = []
        self.DuxDz = []
        self.DuzDx = []
        self.DuzDz = []

        # Scalar fields
        self.rho = []
        self.T = []
        self.p = []
        self.e_t = []

        self.qx = []
        self.qz = []

        # Conserved variables
        self.Q = []

        self.E0 = []
        self.E1 = []
        self.E2 = []
        self.E3 = []

        self.G0 = []
        self.G1 = []
        self.G2 = []
        self.G3 = []

        self.Ev0 = []
        self.Ev1 = []
        self.Ev2 = []
        self.Ev3 = []

        self.Gv0 = []
        self.Gv1 = []
        self.Gv2 = []
        self.Gv3 = []

        self.rhs = []

        self.Q_copy = []

        pass

    def set_arrays(self):

        # Velocity vector components
        self.ux = np.zeros([para.Nx+para.A,para.Nz+1])
        self.uz = np.zeros([para.Nx+para.A,para.Nz+1])

        self.DuxDx = np.zeros([para.Nx+para.A,para.Nz+1])
        self.DuxDz = np.zeros([para.Nx+para.A,para.Nz+1])
        self.DuzDx = np.zeros([para.Nx+para.A,para.Nz+1])
        self.DuzDz = np.zeros([para.Nx+para.A,para.Nz+1])

        # Scalar fields
        self.rho = np.zeros([para.Nx+para.A,para.Nz+1])
        self.T = np.zeros([para.Nx+para.A,para.Nz+1])
        self.p = np.zeros([para.Nx+para.A,para.Nz+1])
        self.e_t = np.zeros([para.Nx+para.A,para.Nz+1])

        self.qx = np.zeros([para.Nx+para.A,para.Nz+1])
        self.qz = np.zeros([para.Nx+para.A,para.Nz+1])

        # Conserved variables
        self.Q = np.zeros([4,para.Nx+para.A,para.Nz+1])

        self.E0 = np.zeros([para.Nx+para.A,para.Nz+1])
        self.E1 = np.zeros([para.Nx+para.A,para.Nz+1])
        self.E2 = np.zeros([para.Nx+para.A,para.Nz+1])
        self.E3 = np.zeros([para.Nx+para.A,para.Nz+1])

        self.G0 = np.zeros([para.Nx+para.A,para.Nz+1])
        self.G1 = np.zeros([para.Nx+para.A,para.Nz+1])
        self.G2 = np.zeros([para.Nx+para.A,para.Nz+1])
        self.G3 = np.zeros([para.Nx+para.A,para.Nz+1])
        
        self.Ev0 = np.zeros([para.Nx+para.A,para.Nz+1])
        self.Ev1 = np.zeros([para.Nx+para.A,para.Nz+1])
        self.Ev2 = np.zeros([para.Nx+para.A,para.Nz+1])
        self.Ev3 = np.zeros([para.Nx+para.A,para.Nz+1])

        self.Gv0 = np.zeros([para.Nx+para.A,para.Nz+1])
        self.Gv1 = np.zeros([para.Nx+para.A,para.Nz+1])
        self.Gv2 = np.zeros([para.Nx+para.A,para.Nz+1])
        self.Gv3 = np.zeros([para.Nx+para.A,para.Nz+1])

        self.rhs = np.zeros([4,para.Nx+para.A,para.Nz+1])

        self.Q_copy = np.zeros([4,para.Nx+para.A,para.Nz+1])

        pass

    def init_hydro(self):

        # with h5py.File("output\\2D_%d.h5" %(199), "r") as f:
        #     # List all groups
        #     # print("Keys: %s" % f.keys())
            
        #     data_T = list(f.keys())[0]
        #     data_p = list(f.keys())[1]
        #     data_ux = list(f.keys())[2]
        #     data_uz = list(f.keys())[3]

        #     # Get the data
        #     self.T = np.array(f[data_T])
        #     self.p = np.array(f[data_p])
        #     self.ux = np.array(f[data_ux])
        #     self.uz = np.array(f[data_uz])
        
        self.T = para.T_l - para.DeltaT*para.Z_mesh

        self.p = (para.C/(para.DeltaT**(para.m+1)))*self.T**(para.m+1)

        self.ux = 0.0001*np.sin(para.X_mesh*np.pi)*np.cos(para.Z_mesh*np.pi*para.A/(self.p/(para.R_gas*self.T)))
        self.uz = - 0.0001*np.cos(para.X_mesh*np.pi)*np.sin(para.Z_mesh*np.pi/(self.p/(para.R_gas*self.T)))

        # self.ux = 0.000001*np.ones([para.Nx+para.A,para.Nz+1])
        # self.uz = np.zeros([para.Nx+para.A,para.Nz+1])

        pass

    def density(self):

        self.rho = self.p/(para.R_gas*self.T)

        pass

    def pressure(self):

        self.p = self.rho*para.R_gas*self.T

        pass

    def specific_energy(self):

        self.e_t = (1/2)*(self.ux**2 + self.uz**2) + para.Cv*self.T + para.g*para.Z_mesh

        pass

    def update_conserved(self):

        self.Q[0] = self.rho
        self.Q[1] = self.rho*self.ux
        self.Q[2] = self.rho*self.uz
        self.Q[3] = self.rho*self.e_t

        pass

    def compute_convective(self):

        self.E0 = self.rho*self.ux
        self.E1 = self.rho*(self.ux**2) + self.p
        self.E2 = self.rho*self.ux*self.uz
        # print('uz1 = ', self.uz)
        self.E3 = (self.rho*self.e_t + self.p)*self.ux

        self.G0 = self.rho*self.uz
        self.G1 = self.rho*self.ux*self.uz
        self.G2 = self.rho*(self.uz**2) + self.p
        self.G3 = (self.rho*self.e_t + self.p)*self.uz
        # print('uz2 = ', self.uz)

        pass

    def compute_viscous(self):

        self.Ev0 = 0
        self.Gv0 = 0

        self.DuxDx = copy(dfx_c(self.ux))
        self.DuxDz = copy(dfz_c(self.ux))

        self.DuzDx = copy(dfx_c(self.uz))
        self.DuzDz = copy(dfz_c(self.uz))

        self.Ev1 = (4/3)*self.DuxDx - (2/3)*self.DuzDz
        self.Ev2 = (self.DuxDz + self.DuzDx)

        self.Gv1 = self.Ev2
        self.Gv2 = (4/3)*self.DuzDz - (2/3)*self.DuxDx

        self.qx = copy(dfx_c(self.T))
        self.qz = copy(dfz_c(self.T))

        self.Ev3 = self.ux*self.Ev1 + \
            self.uz*self.Ev2 + para.K*self.qx
        self.Gv3 = self.ux*self.Gv1 + \
            self.uz*self.Gv2 + para.K*self.qz

        pass

    def viscous_derivative(self):

        self.Ev0 = 0

        self.Ev1 = copy(dfx_c(self.Ev1))
        self.Ev2 = copy(dfx_c(self.Ev2))
        self.Ev3 = copy(dfx_c(self.Ev3))

        self.Gv0 = 0
        self.Gv1 = copy(dfz_c(self.Gv1))
        self.Gv2 = copy(dfz_c(self.Gv2))
        self.Gv3 = copy(dfz_c(self.Gv3))

        pass

    def convective_derivative(self):

        self.E0 = copy(dfx_c(self.E0))
        self.E1 = copy(dfx_c(self.E1))
        self.E2 = copy(dfx_c(self.E2))
        self.E3 = copy(dfx_c(self.E3))

        self.G0 = copy(dfz_c(self.G0))
        self.G1 = copy(dfz_c(self.G1))
        self.G2 = copy(dfz_c(self.G2))
        self.G2 += self.rho*para.g
        self.G3 = copy(dfz_c(self.G3))

        pass

    def compute_rhs(self):

        self.update_conserved()
        self.compute_convective()
        self.convective_derivative()
        self.compute_viscous()
        self.viscous_derivative()

        self.rhs[0] = (self.Ev0 + self.Gv0 - self.E0 - self.G0)
        self.rhs[1] = (self.Ev1 + self.Gv1 - self.E1 - self.G1)
        self.rhs[2] = (self.Ev2 + self.Gv2 - self.E2 - self.G2)
        self.rhs[3] = (self.Ev3 + self.Gv3 - self.E3 - self.G3)

        pass

