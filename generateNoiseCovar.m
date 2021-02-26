function Q = generateNoiseCovar(system_matrix_A, noise_gain_matrix_G,sigma,delta_T)

intfunc = @(tp) expm(system_matrix_A * (delta_T-tp)) * noise_gain_matrix_G ...
    * sigma^2 * noise_gain_matrix_G'* expm(system_matrix_A' * (delta_T-tp));
Q = integral(intfunc,0,delta_T,'ArrayValued',true);
end