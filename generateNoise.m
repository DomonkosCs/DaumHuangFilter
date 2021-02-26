G = [0;1];
A = [0,1;-3*g/2/l,0];
sigma = 1;

dT = 0:0.01:1;
Q = zeros(2,2,numel(dT));
tic
for i = 1:numel(dT)
    Q(:,:,i) = generateNoiseCovar(A,G,sigma,dT(i)); 
end
q22interp = interp1(dT,reshape(Q(2,2,:),[1,numel(dT)]),.8,'spline');
toc

function Q = generateNoiseCovar(system_matrix_A, noise_gain_matrix_G,sigma,delta_T)

intfunc = @(tp) expm(system_matrix_A * (delta_T-tp)) * noise_gain_matrix_G ...
    * sigma^2 * noise_gain_matrix_G'* expm(system_matrix_A' * (delta_T-tp));
Q = integral(intfunc,0,delta_T,'ArrayValued',true);
end