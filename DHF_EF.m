%% 1D test case

% TODO: function descriptions

% plant_function = @(x,k) x/2 + 25*x/(1+ x^2) + 8*cos(1.2*k);
% measurement_function = @(x) x^2/20;
% 
% plant = @(x,k) plant_function(x,k) + normrnd(0,1);
% measurement = @(x) measurement_function(x) + normrnd(0,0.1);

%-------
state_equation = @(xk_prev,k) xk_prev./2 + 25.*xk_prev./(1 + xk_prev.^2) + 8*cos(1.2*k);
output_equation = @(xk) xk.^2 ./ 20;
%-------
%[x,y,z] = pendulumSimulation(0.8,0.01,20,[pi/16,0]);
%plot(x,y)

equations = systemEquationsHandler(state_equation,output_equation);

kInitState = 0;
kDeltaTime = 1;
kProcessNoiseVarianceQ = 10;
kMeasurementNoiseVarianceR = 1;
kParticleCount = 130;
kMaxProcessTime = 1000;
kMonteCarloRuns = 1000;

% [state_array,measurement_array] = simulateProcess(equations,kInitState,kDeltaTime,kMaxProcessTime, ...
%             kProcessNoiseVarianceQ,kMeasurementNoiseVarianceR);
% [dhf_filtered_state,ekf_filtered_state] = exactFlowFilter(equations,kParticleCount,kMaxProcessTime,kDeltaTime,...
%             kProcessNoiseVarianceQ,kMeasurementNoiseVarianceR,measurement_array);
% plot(dhf_filtered_state)
% hold on
% plot(state_array)
% plot(ekf_filtered_state)
% legend('dhf','true','ekf')

% filter_RMSE = zeros(kMonteCarloRuns,kMaxProcessTime);

for i=200
    parfor j = 1:kMonteCarloRuns
        [state_array,measurement_array] = ...
            simulateProcess(equations,kInitState,kDeltaTime,i, ...
            kProcessNoiseVarianceQ,kMeasurementNoiseVarianceR);
        [dhf_filtered_state,ekf_filtered_state] = exactFlowFilter(equations,kParticleCount,i,kDeltaTime,...
            kProcessNoiseVarianceQ,kMeasurementNoiseVarianceR,measurement_array);
        filter_RMSE(j,i) = sqrt(mean((state_array-ekf_filtered_state).^2));
    end
end
filter_RMSE_mean = mean(filter_RMSE);
%figure(2)
%plot(filter_RMSE_mean);


function systemEquations = systemEquationsHandler(state_equation,output_equation)
    syms xk xk_prev k
    f(xk_prev,k) = sym(state_equation);
    h(xk) = sym(output_equation);
    
    white_noise_func = @(X,array_size) normrnd(0,X,array_size);
    plant_with_noise_func = @(xk_prev,k,Q) state_equation(xk_prev,k) ...
        + white_noise_func(sqrt(Q),size(xk_prev));
    sensor_with_noise_func = @(xk,R) output_equation(xk) ...
        + white_noise_func(sqrt(R),size(xk));

    linear_plant = matlabFunction(jacobian(f,xk_prev));
    linear_sensor = matlabFunction(jacobian(h,xk));
    systemEquations.plant_with_noise = plant_with_noise_func;
    systemEquations.plant_noiseless = state_equation;
    systemEquations.sensor_with_noise = sensor_with_noise_func;
    systemEquations.sensor_noiseless = output_equation;
    systemEquations.linear_plant = linear_plant;
    systemEquations.linear_sensor = linear_sensor;
end

function [state_array,measurement_array] = simulateProcess(system_equations, ...
        init_state,delta_time,end_time, ...
        process_noise_variance_Q, measurement_noise_variance_R)
    
    state_array = generateState(system_equations,init_state, ...
        delta_time,end_time,process_noise_variance_Q);
    measurement_array = system_equations.sensor_with_noise(state_array, ...
        measurement_noise_variance_R);

end

function state_x = generateState(system_equations,init_state,delta_time, ...
        end_time,process_noise_variance_Q)
    
    kStartTime = 0;
    time_eval_points = kStartTime:delta_time:end_time;
    state_x = zeros(1,numel(time_eval_points));
    
    state_x(1) = init_state;
    for i = 1:end_time/delta_time
        current_time = (i - 1) * delta_time;
        current_state = state_x(i);
        next_time = current_time + delta_time;
        state_x(i+1) = ... % next_state
            system_equations.plant_with_noise(current_state,next_time,process_noise_variance_Q);
    end
end

function [EKF_mean_m,filtered_state] = exactFlowFilter(system_equations,particle_count,end_time,delta_time, ...
        process_noise_variance_Q,measurement_noise_variance_R,measurement_z)
    
    kLambdaDivPoints = 11;
    kInitialDistributionCovar = 1;
    
    d_lambda = 1/(kLambdaDivPoints-1);
    time_div_points = end_time/delta_time+1;

    filtered_state = zeros(1,time_div_points);
    
    % rows: particles; coloumns: division in pseudo time; layers: division in time
    state_array = zeros(particle_count,kLambdaDivPoints,time_div_points);
    
    predicted_EKF_mean_m = zeros(1,time_div_points);
    predicted_EKF_covar_P = zeros(1,time_div_points);
    EKF_mean_m = zeros(1,time_div_points);
    EKF_covar_P = zeros(1,time_div_points);
    
    % initialization from a gaussian, centered around the initial cond (if know)
    EKF_covar_P(1) = kInitialDistributionCovar;
    particle_vector_size = [particle_count,1];
    state_array(:,1,1) = normrnd(0,sqrt(EKF_covar_P(1)),particle_vector_size);
    
    filtered_state(1) = mean(state_array(:,1,1));
    EKF_mean_m(1) = filtered_state(1);
    
    % iterate through time
    for i = 2:time_div_points
        state_array(:,1,i) = system_equations.plant_with_noise(state_array(:,1,i-1),i, ...
            process_noise_variance_Q);
        filtered_state(i) = mean(state_array(:,1,i));
        [predicted_EKF_mean_m(i),predicted_EKF_covar_P(i)] = ...
            EKFPrediction(system_equations,EKF_mean_m(i-1),EKF_covar_P(i-1),process_noise_variance_Q,i);
        % iterate through the pseudo time lambda
        for j = 1:kLambdaDivPoints
            lambda = j * d_lambda;
            % iterate through the particles
            for k = 1:particle_count
                dx_dlambda = calculateSlope(system_equations,lambda,state_array(k,j,i),filtered_state(i), ...
                    measurement_noise_variance_R,predicted_EKF_covar_P(i),measurement_z(i));
                state_array(k,j+1,i) = state_array(k,j,i) + dx_dlambda * d_lambda; 
            end
            filtered_state(i) = mean(state_array(:,j+1,i)); % re evaluate the avg
        end
        [EKF_mean_m(i),EKF_covar_P(i)] = EKFUpdate(system_equations,predicted_EKF_mean_m(i), ...
            predicted_EKF_covar_P(i),measurement_noise_variance_R,measurement_z(i));
    end
end % TODO: get rid of the for loops -> vectorize

function dx_dlambda = calculateSlope(system_equations,lambda,currentState,state_avg, ...
        measurement_noise_variance_R,predicted_EKF_covar_P,measurement_z)
    % Calculate the linearized measurement matrix for each particle
    % dx/dlambda = B*x + b
    
    linear_measurement_matrix_H = system_equations.sensor_noiseless(currentState); % TODO function
    B = -1/2*predicted_EKF_covar_P*linear_measurement_matrix_H' ...
        /(lambda*linear_measurement_matrix_H*predicted_EKF_covar_P ...
        *linear_measurement_matrix_H' + measurement_noise_variance_R) * linear_measurement_matrix_H;
    b = (eye(size(B)) + 2*lambda*B) ...
        * ((eye(size(B)) + lambda*B)*predicted_EKF_covar_P*linear_measurement_matrix_H'...
        / measurement_noise_variance_R*measurement_z + B*state_avg);
    
    dx_dlambda = B*currentState + b;
end

function [predicted_EKF_mean_m,predicted_EKF_covar_P] = ...
        EKFPrediction(system_equations,EKF_mean_m, EKF_covar_P, ...
                      process_noise_variance_Q,time)
    
    predicted_EKF_mean_m = system_equations.plant_noiseless(EKF_mean_m,time);
    linear_plant_eval = system_equations.linear_plant(EKF_mean_m);
    predicted_EKF_covar_P = linear_plant_eval*EKF_covar_P*linear_plant_eval' ...
        + process_noise_variance_Q; 
    % TODO: make sure that the linearized system matrix is evaluated at the
    % EKF_mean_m point.
    % TODO: sometimes there is matrix multiplication in the process noise
end

function [EKF_mean_m,EKF_covar_P] = EKFUpdate(system_equations,predicted_EKF_mean_m, ... 
        predicted_EKF_covar_P,measurement_noise_variance_R,measurement_z)
    
    linear_sensor_eval = system_equations.linear_sensor(predicted_EKF_mean_m);
    kalman_gain_K = predicted_EKF_covar_P * linear_sensor_eval' ...
        /(linear_sensor_eval*predicted_EKF_covar_P*linear_sensor_eval' ...
        + measurement_noise_variance_R);
    EKF_mean_m = predicted_EKF_mean_m + kalman_gain_K ...
        *(measurement_z - system_equations.sensor_noiseless(predicted_EKF_mean_m));
    EKF_covar_P = (1-kalman_gain_K*linear_sensor_eval)*predicted_EKF_covar_P;
end

