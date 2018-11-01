# Waddup kids!
import numpy as np
import math
import time as tm
import matplotlib.pyplot as plt

N = 6  # at the moment only even values of N allowed due to initialisation
w = 4
dt = 0.1
dx = 0.2
rho = 1.05
v = 2.05
P = 3.05
gamma = 1.4


def make_mesh():  # Function that makes the N-Dimensional mesh filled with zero values.
    # first row  : we work in this "point"
    # second row : value of rho here
    # third row  : value of v here
    # fourth row : value of P here
    mesh = np.zeros([4, N + 2])
    return mesh


def initialisation(array, rho_1, v_1, P_1, rho_2, v_2, P_2):
    for i in range(N + 2):
        if i < (N + 2) / 2:
            array[0, i] = i
            array[1, i] = rho_1
            array[2, i] = v_1
            array[3, i] = P_1
        else:
            array[0, i] = i
            array[1, i] = rho_2
            array[2, i] = v_2
            array[3, i] = P_2
    return array


def comp_diagonal(initialised_mesh, i, gamma):  # we also need to calculate the diagonal on each timestep
    rho = initialised_mesh[1, i]
    v = initialised_mesh[2, i]
    P = initialised_mesh[3, i]

    diagonal = np.array( [[v, 0, 0],
                         [0, v + math.sqrt((gamma * P) / rho), 0],
                         [0, 0, v - math.sqrt( (gamma * P) / rho)]])
    return diagonal

def comp_A(initialised_mesh, i, gamma):
    rho = initialised_mesh[1, i]
    v = initialised_mesh[2, i]
    P = initialised_mesh[3, i]
    A = np.array([[1, rho / (gamma * P), rho / (gamma * P)],
                  [0, math.sqrt(1 / (gamma * P * rho)), - math.sqrt(1 / (gamma * P * rho))],
                  [0, 1, 1]])
    return A

def comp_inverse(A): #gives the A^-1 matrix needed to solve everything
    A_inverse = np.linalg.inv(A)
    return A_inverse

def create_K(initialised_mesh,i):
    K = np.zeros([3,1])
    K[0,0] = initialised_mesh[1, i]
    K[1,0] = initialised_mesh[2, i]
    K[2,0] = initialised_mesh[3, i]
    return K

def solver(mesh,dt,dx, rho, v , P , gamma ): #we solve for one time step !
   #for j in range(0,1,dt): #go over the timesteps
       for i in range(1, N+1): #go over all points in the mesh
           D = comp_diagonal(mesh,i, gamma)
           A_1 =comp_A(mesh,i,gamma)
           A = comp_inverse(A_1)
           K = create_K(mesh, i)
           W = np.matmul(A , K)
           u = {}
           for k in range(3):
               u['u_%01d_%01d'%(i,k)] =  mesh[k, i] - D[k,k] * dt/dx * (mesh[k, i]-mesh[k, i-1])
           print(u)
          # for k in range(3)


           #u_i_j_dt =j


"""def initialization (x, w):
def computetimestep (w, dt):
def getcmax (w, cmax):
def integrateintime(x, w, dt):
def getbc (x, w):
def get_flux (w, lam, f_upwind):
def cons_to_diag (w, wdiag):
def get_kinv (v, csound, Kinv):
def diag_to_cons (wdiag, w):
def get_k (v, csound, K):
def get_eigenval (w, lam):
def primitive (w, wprim):
def cv_and_conquer (x, w):
def chrono (start, start_mess):
def followup (message):"""

A = make_mesh()
initialisation(A, rho, v, P, 4, 5, 6)
#print (create_K(A, 1))
# print (compute_inverse(1,1,1,1))
solver(A,dt,dx, rho, v , P , gamma )

"""for i in range(1,N+2):
    print(i)"""