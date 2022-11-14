import numpy as np
from pathlib import Path
import para
import h5py
from derivative import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from pylab import rcParams

path=Path.cwd()

# Energy plotting in 2D

# with h5py.File("output\E&Nu.h5", "r") as k:

#     print("Keys: %s" % k.keys())

#     data_E = list(k.keys())[0]
#     data_K = list(k.keys())[1]
#     data_e = list(k.keys())[2]
#     data_Nu = list(k.keys())[3]

#     E = np.array(k[data_E])
#     K = np.array(k[data_K])
#     e = np.array(k[data_e])
#     Nu = np.array(k[data_Nu])

#     print(len(E),len(e),len(para.t))

#     plt.plot(para.t[:-1],Nu)
#     # plt.plot(para.t[:-1],E,label='Total energy')
#     # plt.plot(para.t[:-1],K,label='Total kinetic energy')
#     # plt.plot(para.t[:-1],e,label='Total internal energy')

#     plt.xlabel(r'$t$')
#     # plt.ylabel(r'$Energy$')
#     plt.ylabel(r'$Nu$')
#     plt.xlim(0,3.5)
#     # plt.ylim(0,250)
#     plt.grid(color = 'k', linestyle = '--', linewidth = 0.5)
#     plt.legend()
#     # plt.title('Energy Plot')
#     plt.title('Nusselt number vs $t$ plot')
#     plt.tight_layout()
#     plt.savefig('postprocess\\Nu_vs_t.png')
#     # plt.savefig('postprocess\\E_vs_t.png')

# Fields plotting in 2D

# hra = np.array([0,72])
# for i in hra:

#     with h5py.File("output\2D_%d.h5" %(i), "r") as f:
#         # List all groups
#         # print("Keys: %s" % f.keys())
        
#         data_T = list(f.keys())[0]
#         data_p = list(f.keys())[1]
#         data_ux = list(f.keys())[2]
#         data_uz = list(f.keys())[3]

#         T = np.array(f[data_T])
#         p = np.array(f[data_p])
#         rho = p/(para.R_gas*T)
#         ux = np.array(f[data_ux])
#         uz = np.array(f[data_uz])

#         dTdz = dfz_c(T)
#         dTdz_h = np.sum(dTdz,axis=1)/(para.Nx+1)
#         T_ha = np.sum(dTdz_h)/(para.Nz+1)
#         F_a = para.g*para.K/para.Cp
#         F_c = para.K*(para.epsilon - 1)/para.Lz
#         F_t = para.K*(para.epsilon + T_ha)
#         N = (F_t - F_a)/(F_c - F_a)

#         print ('Nusselt number at bottom plate at %f = ' %(np.round(i,2)), N)

#         T_h = np.sum(T,axis = 0)/(para.Nx+1)
#         p_h = np.sum(p,axis = 0)/(para.Nx+1)
#         rho_h = np.sum(rho,axis = 0)/(para.Nx+1)
        
#         plt.plot(para.Z,T_h,label='T at t = %f' %(np.round(para.tstep[i],2)))
#         # plt.plot(para.Z,para.T_l - para.epsilon*para.Z,label='T_a')
#         # plt.plot(para.Z,p_h/p[:,0],label='p')
#         # plt.plot(para.Z,rho_h/rho[:,0],label='rho')
#         plt.xlabel(r'$Z$')
#         # plt.ylabel(r'$$')
#         plt.xlim(0,1)
#         # plt.ylim(0,9000)
#         plt.grid(color = 'k', linestyle = '--', linewidth = 0.5)
#         plt.legend()
#         plt.title('Horizontal averages')
#         plt.tight_layout()
#         plt.savefig('hori_avg.png')
#         # plt.show()

#         pass
    
#     pass

i = 15
while i < para.Nn:

    with h5py.File("output\\2D_%d.h5" %(i), "r") as f:
        # List all groups
        # print("Keys: %s" % f.keys())
        
        data_T = list(f.keys())[0]
        data_p = list(f.keys())[1]
        data_ux = list(f.keys())[2]
        data_uz = list(f.keys())[3]


        # Get the data
        T = np.array(f[data_T])
        p = np.array(f[data_p])
        rho = p/(para.R_gas*T)
        ux = np.array(f[data_ux])
        uz = np.array(f[data_uz])
        e_t = (1/2)*(ux**2 + uz**2) + para.Cv*T

        # Vector Plot for velocity vector in 2D box
        fig2 = plt.figure(figsize=(16,8))
        ax2 = fig2.add_subplot(aspect = 'equal')
        # c2 = ax2.pcolormesh(para.X_mesh,para.Z_mesh,T)
        c2 = ax2.contourf(para.X_mesh,para.Z_mesh,T,100,cmap=cm.coolwarm)
        ax2.set_aspect(aspect = 1)
        divider = make_axes_locatable(ax2)
        cax2 = divider.append_axes('right', size='5%', pad=0.05)
        fig2.colorbar(c2, cax=cax2,label='$T$')
        ax2.quiver(para.X_mesh[::6,::6],para.Z_mesh[::6,::6],ux[::6,::6],uz[::6,::6],color = 'k',label='$\vec{u}$',linewidths = 1)
        ax2.set_xlim(0,para.A)
        ax2.set_ylim(0,1)
        ax2.set_xlabel('$X$')
        ax2.set_ylabel('$Z$')
        ax2.set_title('t = %f' %(para.tstep[i]))

    plt.savefig('postprocess\\%d.png' %(i))
    # plt.show()

    i = i+1

    pass

