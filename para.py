from re import T
import numpy as np


# ----------------------------- Time variables ----------------------------- #

tinit = 0                                           # Initial time
tfinal = 10                                      # Final time
dt = 1e-6                                           # Single time step
t = np.arange(tinit,tfinal+dt,dt)                   # Time axis

Nn = 500
tstep = np.linspace(tinit,tfinal,Nn+1)

n1=1                                                # Used for rounding time in time-stepping
n2=10 
while dt*n2!=1:
    n1=n1+1
    n2=n2*10

# -------------------------------------------------------------------------- #


# ----------------------- Characteristics parameters ----------------------- #

alpha = 1.5                                 # Adiabatic index
gamma = (1+alpha)/alpha                     # Ratio of specific heat
A = 2                                       # Aspect ratio
Pr = 1                                      # Prandtl number
DeltaT = 1                                  # Normalised layer thickness parameter
m = 1.4                                     # Polytropic index
Ra = 1e3                                    # Rayleigh number

# -------------------------------------------------------------------------- #


# ----------------------------- Grid variables ----------------------------- #

Nz = 128                                    # Number of grid points in z-direction
Nx = A*Nz                                   # Number of grid points in x-direction

Lz = 1                                      # Length of box in z-direction
Lx = A                                      # Length of box in x-direction

dz = Lz/Nz                                  # Length between two consecutive grid points in z-direction         
dx = Lx/Nx                                  # Length between two consecutive grid points in x-direction

Z = np.arange(0,(Nz+1))*dz                  # Array consisting of points in z-direction
X = np.arange(0,(Nx+A))*dx                  # Array consisting of points in x-direction

X_mesh, Z_mesh = np.meshgrid(X, Z,indexing = 'ij')       # Meshgrids

# -------------------------------------------------------------------------- #


# --------------------------- Related parameters --------------------------- #

g = (1+alpha)*Ra/(Pr*DeltaT*(alpha-m))      # Acceleration due to gravity
R_gas = g/(DeltaT*(m+1))                    # Real gas constant
K = (1+alpha)*R_gas/Pr                      # Thermal conductivity
T_l = DeltaT + 1                            # Fixed temperature at lower plate
C = R_gas/DeltaT**(m+1)                     # Mass constant
Cp = (1+alpha)*R_gas                        # Specific heat capacity at constant pressure
Cv = alpha*R_gas                            # Specific heat capacity at constant volume

Chi = (C/(R_gas*DeltaT**(m+1)))*T_l**m      # Initial density constrast (rho_l/rho_u) 

v_sound = np.sqrt(gamma*R_gas*T_l)          # Speed of sound near bottom plate

t_cfl = dx/v_sound                          # CFL conditon time

if dt > t_cfl:

    print ('Error! dt is larger than CFL time, please lower dt.')

    exit()

# -------------------------------------------------------------------------- #


print ('g = ', g, 'R_gas = ', R_gas, 'K = ', K, 'C = ', C, 'Chi = ', Chi, ' Cp = ', Cp, 'v_sound = ', v_sound, 't_cfl = ', t_cfl)