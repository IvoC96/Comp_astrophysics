syms rho v gamma P
A = [v rho 0; 0 v 1/rho; 0 gamma*P v];
e = eig(A)
[D, E] = eig(A);
F= simplify(D) 
G = simplify(E)
