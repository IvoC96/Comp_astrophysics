syms rho v gamma P
% A = [v rho 0; 0 v 1/rho; 0 gamma*P v];
% e = eig(A)
% [D, E] = eig(A);
% F= simplify(D) 
% G = simplify(E)

H = [1 rho./(gamma*P) rho./(gamma*P); 0 sqrt(1/(gamma* rho* P)) -sqrt(1/(gamma* rho* P)); 0 1 1 ];
I = inv(H)
