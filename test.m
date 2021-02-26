syms delta_T tp sigma positive real
 
%generateNoiseCovar([0,1;0,0],[0;1],sigma,dt)

system_matrix_A = [0,1,0;0,0,1;0,0,0];
A = [0,1;-3*g/2/l,0];
noise_gain_matrix_G = [0;0;1];


int(expm(system_matrix_A * (delta_T-tp)) * noise_gain_matrix_G ...
    * sigma^2 * noise_gain_matrix_G'* expm(system_matrix_A' * (delta_T-tp)),tp,0,delta_T)