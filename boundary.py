import para
from compressible import Compressible

def imposeBC_u(compress = Compressible()):
    # Free-slip boundary condition for velocity

    compress.ux[0,:] = compress.ux[-1,:] = 0 

    compress.ux[:,0] = compress.ux[:,1] 
    compress.ux[:,-1] = compress.ux[:,-2]  

    compress.uz[0,:] = compress.uz[1,:]
    compress.uz[-1,:]  = compress.uz[-2,:]

    compress.uz[:,0] = compress.uz[:,-1] = 0

    pass

def imposeBC_T(compress = Compressible()):
    # Fixed boundary on vertical walls, adiabatic on horizontal walls for temperature

    compress.T[0,:] = compress.T[1,:] 
    compress.T[-1,:] = compress.T[-2,:]

    compress.T[:,0] = para.DeltaT + 1
    compress.T[:,-1] = 1

    pass

def imposeBC_p(compress = Compressible()):
    # Neumann boundary condition for pressure 

    compress.p[:,0] = compress.p[:,1]
    compress.p[:,-1] = compress.p[:,-2]

    compress.p[0,:] = compress.p[1,:]
    compress.p[-1,:] = compress.p[-2,:] 

    pass

def boundary(compress = Compressible()):

    imposeBC_u(compress)
    imposeBC_T(compress)
    imposeBC_p(compress)

    pass