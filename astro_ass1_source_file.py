# Waddup kids!
import numpy as np
import math
import time as tm
import matplotlib.pyplot as plt

N = 6 #at the moment only even values of N allowed due to initialisation
w = 4
dt = 5
rho = 1
v = 2
P = 3
gamma = 1.4


def make_mesh():  # Function that makes the N-Dimensional mesh filled with zero values.
    # first row  : we work in this "point"
    # second row : value of rho here
    # third row  : value of v here
    # fourth row : value of P here
    mesh = np.zeros([4, N + 2])
    return mesh

def initialisation(array, rho_1 , v_1 , P_1,rho_2, v_2, P_2):
    for i in range(N+2):
        if i < (N+2)/2:
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

def compute_diagonal(rho, v, P, gamma):  # we also need to calculate the diagonal on each timestep
    diagonal = np.array([v, 0, 0], [[0, v + math.sqrt((gamma * P) / rho), 0], [0, 0, v - math.sqrt((gamma * P) / rho)]])
    return diagonal


def compute_inverse(rho, v, P, gamma):
    A = np.array([[1, rho / (gamma * P), rho / (gamma * P)],
                 [0, math.sqrt(1 / (gamma * P * rho)), - math.sqrt(1 / (gamma * P * rho))],
                 [0, 1, 1]])
    A_inverse = np.linalg.inv(A)
    return A_inverse


"""def initialization (x, w): #each point in the mesh should be initiated with the values as given by the user.
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
print (initialisation(A, rho , v , P,4, 5, 6))

#print (compute_inverse(1,1,1,1))
