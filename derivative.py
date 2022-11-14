import numpy as np 
import para
from numpy import copy

# Central-difference method

def dfx_c(f):
    # Derivative of array 'f' w.r.t. x

    a = ((np.roll(f,-1,axis=0) - np.roll(f,1,axis=0))/(2*para.dx))

    a[0,:] = (f[1,:] - f[0,:])/para.dx
    a[-1,:] = (f[-1,:] - f[-2,:])/para.dx

    return a

def dfz_c(f):
    # Derivative of array 'f' w.r.t. z

    a = ((np.roll(f,-1,axis=1) - np.roll(f,1,axis=1))/(2*para.dz))

    a[:,0] = (f[:,1] - f[:,0])/para.dz
    a[:,-1] = (f[:,-1] - f[:,-2])/para.dz

    return a

# 6th order Compact Scheme

def thomas(d):

    N = len(d)

    a = (1/3)*np.ones(N,dtype='float64')
    a[0] = 0.0
    a[1] = 2/11
    a[-2] = 2/11
    a[-1] = 5.0
    c = (1/3)*np.ones(N,dtype='float64')
    c[0] = 5.0
    c[1] = 2/11
    c[-2] = 2/11
    c[-1] = 0.0

    # print (a,b,c)

    cp = np.zeros(N,dtype='float64') # store tranformed c or c'
    dp = np.zeros(N,dtype='float64') # store transformed d or d'
    X = np.zeros(N,dtype='float64') # store unknown coefficients
    
    # Perform Forward Sweep
    # Equation 1 indexed as 0 in python
    cp[0] = c[0] 
    dp[0] = d[0]
    # Equation 2, ..., N (indexed 1 - N-1 in Python)
    for i in np.arange(1,(N),1):
        dnum = 1 - a[i]*cp[i-1]
        cp[i] = c[i]/dnum
        dp[i] = (d[i]-a[i]*dp[i-1])/dnum
    
    # Perform Back Substitution
    X[(N-1)] = dp[N-1]  # Obtain last xn 

    for i in np.arange((N-2),-1,-1):  # use x[i+1] to obtain x[i]
        X[i] = (dp[i]) - (cp[i])*(X[i+1])
    
    d[:] = X

    pass

def x_derv(B):

    f = np.ones([para.Nx+1,para.Nz+1])

    f[0,:] = (B[1,:] - B[0,:])/para.dz
    f[-1,:] = (B[-1,:] - B[-2,:])/para.dz

    A = B[1:-1,:]

    temp = (1/9)*(np.roll(A,-2,axis=0) - np.roll(A,2,axis=0))/(4*para.dx) + \
            (14/9)*(np.roll(A,-1,axis=0) - np.roll(A,1,axis=0))/(2*para.dx)

    temp[0,:] = ((-197/60)*A[0,:] + (-5/12)*A[1,:] + (5)*A[2,:] + (-5/3)*A[3,:] \
                     + (5/12)*A[4,:]+ (-1/20)*A[5,:])/para.dx
    temp[1,:] = ((-20/33)*A[0,:] + (-35/132)*A[1,:] + (34/33)*A[2,:] + (-7/33)*A[3,:] \
                     + (2/33)*A[4,:]+ (-1/132)*A[5,:])/para.dx 

    temp[-2,:] = ((20/33)*A[-1,:] + (35/132)*A[-2,:] + (-34/33)*A[-3,:] + (7/33)*A[-4,:] \
                     + (-2/33)*A[-5,:]+ (1/132)*A[-6,:])/para.dx
    temp[-1,:] = ((197/60)*A[-1,:] + (5/12)*A[-2,:] + (-5)*A[-3,:] + (5/3)*A[-4,:] \
                     + (-5/12)*A[-5,:]+ (1/20)*A[-6,:])/para.dx

    for i in range(0,para.Nz+1):
        thomas(temp[:,i])
        pass

    f[1:-1,:] = temp

    return f

def z_derv(B):

    f = np.ones([para.Nx+1,para.Nz+1])

    f[:,0] = (B[:,1] - B[:,0])/para.dz
    f[:,-1] = (B[:,-1] - B[:,-2])/para.dz

    A = B[:,1:-1]

    temp = (1/9)*(np.roll(A,-2,axis=1) - np.roll(A,2,axis=1))/(4*para.dz) + \
            (14/9)*(np.roll(A,-1,axis=1) - np.roll(A,1,axis=1))/(2*para.dz)

    temp[:,0] = ((-197/60)*A[:,0] + (-5/12)*A[:,1] + (5)*A[:,2] + (-5/3)*A[:,3] \
                     + (5/12)*A[:,4]+ (-1/20)*A[:,5])/para.dz
    temp[:,1] = ((-20/33)*A[:,0] + (-35/132)*A[:,1] + (34/33)*A[:,2] + (-7/33)*A[:,3] \
                     + (2/33)*A[:,4]+ (-1/132)*A[:,5])/para.dz 

    temp[:,-2] = ((20/33)*A[:,-1] + (35/132)*A[:,-2] + (-34/33)*A[:,-3] + (7/33)*A[:,-4] \
                     + (-2/33)*A[:,-5]+ (1/132)*A[:,-6])/para.dz
    temp[:,-1] = ((197/60)*A[:,-1] + (5/12)*A[:,-2] + (-5)*A[:,-3] + (5/3)*A[:,-4] \
                     + (-5/12)*A[:,-5]+ (1/20)*A[:,-6])/para.dz

    for i in range(0,para.Nx+1):
        thomas(temp[i,:])
        pass
    
    f[:,1:-1] = temp

    return f