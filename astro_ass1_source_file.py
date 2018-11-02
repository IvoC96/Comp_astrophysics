# Waddup kids!
import numpy as np
import math
import time as tm
import matplotlib.pyplot as plt
x_min = -0.5
x_max = 0.5
N = 20


dt = 0.1
dx = (x_max - x_min)/(N+1)

gamma = 1.4
rho_1 = 8
v_1 = 0
P_1 = 8/gamma

rho_2 = 1
v_2 = 0
P_2 = 1


def make_mesh():
    # Function that makes the N-Dimensional mesh filled with zero values and for
    # the moment being also 2 ghostcells (1 on each side)

    # first row  : we work in this "point"
    # second row : value of rho here
    # third row  : value of v here
    # fourth row : value of P here
    mesh = np.zeros([4, N + 2])
    return mesh


def initial_initialisation(array, rho_1, v_1, P_1, rho_2, v_2, P_2):
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


def initialisation(array, rho, v, P, i):
    array[1, i] = rho[i]
    array[2, i] = v[i]
    array[3, i] = P[i]
    return array


def comp_diagonal(mesh, i, gamma):  # we also need to calculate the diagonal on each timestep
    rho = mesh[1, i]
    v = mesh[2, i]
    P = mesh[3, i]
    diagonal = np.array([[v, 0, 0],
                         [0, v + math.sqrt((gamma * P) / rho), 0],
                         [0, 0, v - math.sqrt((gamma * P) / rho)]])
    """print("in comp_diagonal de diagonal is: ")
    print(diagonal)"""
    return diagonal


def comp_A(mesh, i, gamma):
    rho = mesh[1, i]
    v = mesh[2, i]
    P = mesh[3, i]
    A = np.array([[1, rho / (gamma * P), rho / (gamma * P)],
                  [0, -math.sqrt(1 / (gamma * P * rho)), math.sqrt(1 / (gamma * P * rho))],
                  [0, 1, 1]])
    return A


def comp_inverse(A):  # gives the A^-1 matrix needed to create the vector
    A_inverse = np.linalg.inv(A)
    return A_inverse


def create_K(mesh, i):
    K = np.zeros([3, 1])
    K[0, 0] = mesh[1, i]
    K[1, 0] = mesh[2, i]
    K[2, 0] = mesh[3, i]
    return K


def solver(mesh, dt, dx):  # we solve for one time step !
    #print(mesh)
    U = np.zeros([3, 1])
    rhos = []
    vs = []
    Ps = []
    gamma = 1.4
    for j in np.arange(0.0, 1.0, 1):  # go over the timesteps
        for i in range(1, N + 1):  # go over all points in the mesh
            D = comp_diagonal(mesh, i, gamma)
            A = comp_A(mesh, i, gamma)
            A_inverse= comp_inverse(A)
            K = create_K(mesh, i)
            """print("in solver wordt de diagonal:")
            print(D)"""

            W = np.matmul(A_inverse, K)

            D_2 = comp_diagonal(mesh, i - 1, gamma)
            A_2 = comp_A(mesh, i - 1, gamma)
            A_inverse_2 = comp_inverse(A_2)
            K_2 = create_K(mesh, i - 1)
            W_2 = np.matmul(A_inverse_2, K_2)

            D_3 = comp_diagonal(mesh, i + 1, gamma)
            A_3 = comp_A(mesh, i + 1, gamma)
            A_inverse_3 = comp_inverse(A_3)
            K_3 = create_K(mesh, i + 1)
            W_3 = np.matmul(A_inverse_3, K_3)
            print(dt/dx)
            for k in range(3):
                U[k, 0] = W[k, 0] - (dt/dx) * (
                            (max(D[k, k], 0) * (W[k, 0] - W_2[k, 0])) + (min(D[k, k], 0) * (W_3[k, 0] - W[k, 0])))
            # beforehand, U was equal to A^-1 * K , but we want K for the pure values

            U = np.matmul(A, U)
            rho = U[0, 0]
            v = U[1, 0]
            P = U[2, 0]
            rhos.append(rho)
            vs.append(v)
            Ps.append(P)

        # important we reinitalise the mesh after each time step

        for m in range(N):
            if m == 0:
                mesh[1, m] = rhos[m + 1]
                mesh[2, m] = vs[m + 1]
                mesh[3, m] = Ps[m + 1]
            elif m == N + 1:
                mesh[1, m] = rhos[m - 1]
                mesh[2, m] = vs[m - 1]
                mesh[3, m] = Ps[m - 1]
            else:
                mesh = initialisation(mesh, rhos, vs, Ps, m)
        print(mesh)

        # U = 1D - vector with U[1] = rho, U[2] = v and U[3] = P


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
B = initial_initialisation(A, rho_1, v_1, P_1, rho_2, v_2, P_2)
# print (create_K(A, 1))
# print (compute_inverse(1,1,1,1))

solver(B, dt, dx)
