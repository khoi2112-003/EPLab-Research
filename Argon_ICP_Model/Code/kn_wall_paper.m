function k_wall = kn_wall_paper(Tg, geom, n_g, gamma_X, sigma_nn)
% Paper Eq. (8)

kB  = 1.380649e-23;  MAr = 39.948*1.66053906660e-27;

R = geom.R;  L = geom.L;
A = 2*pi*R^2 + 2*pi*R*L;     
V = pi*R^2*L;                

% neutral thermal speed & mean free path
v_th    = sqrt(8*kB*Tg/(pi*MAr));                    
lambda_n = 1 / max(n_g*max(sigma_nn,1e-25),1e-6);  

% diffusion coefficient 
D_n = (1/3) * v_th * lambda_n;                        

% effective diffusion length 
Lambda_n = 1 / sqrt( (pi/L)^2 + (2.405/R)^2 );        

v_g = v_th;                                           

% assemble k_wall
term_diff = (Lambda_n^2) / max(D_n,1e-30);            
term_surf = (2*V*(2 - gamma_X)) / max(A*v_g*gamma_X,1e-30); 
k_wall    = 1 / (term_diff + term_surf);              
end