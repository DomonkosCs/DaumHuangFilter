%% 1D test case

% plant_function = @(x,k) x/2 + 25*x/(1+ x^2) + 8*cos(1.2*k);
% measurement_function = @(x) x^2/20;
% 
% plant = @(x,k) plant_function(x,k) + normrnd(0,1);
% measurement = @(x) measurement_function(x) + normrnd(0,0.1);

syms f(xk,k) h(xk)
%-------
f(xk,k) = xk/2 + 25*xk/(1+ xk^2) + 8*cos(1.2*k);
h(xk) = xk^2/20;
%-------

sensor_with_noise = h(xk) + normrnd(0,0.1);
plant_with_noise = f(xk,k) + normrnd(0,1);

linearized_plant = matlabFunction(jacobian(f,xk));
linearizer_sensor = matlabFunction(jacobian(h,xk));

testFcn(10,10,linearized_plant)
function test = testFcn(x,k,fcn)
    test = fcn(x,k);
end

function process
kInitState = 0;
kStartTime = 0;
kDeltaTime = 1;
kProcessNoiseVarianceQ = 1;
kMeasurementNoiseVarianceR = 0.1;
kParticleCount = 30;
kMaxProcessTime = 30;
kMonteCarloRuns = 100;

filter_RMSE = zeros(kMonteCarloRuns,kMaxProcessTime);

for i=1:kMaxProcessTime
    parfor j = 1:kMonteCarloRuns
        [state_array,measurement_array] = ...
            simulateProcess(kInitState,kStartTime,kDeltaTime,i, ...
            kProcessNoiseVarianceQ,kMeasurementNoiseVarianceR);
        filtered_state = exactFlowFilter(kParticleCount,i,kDeltaTime,...
            kProcessNoiseVarianceQ,kMeasurementNoiseVarianceR,measurement_array);
        filter_RMSE(j,i) = sqrt(mean((state_array-filtered_state).^2));
        % TODO plot
    end
end
filter_RMSE_mean = mean(filter_RMSE);
plot(filter_RMSE_mean);
end
function [state_array,measurement_array] = simulateProcess(init_state, ...
        start_time,delta_time,end_time, ...
        process_noise_variance_Q, measurement_noise_variance_R)
    
    state_array = generateState(init_state,start_time,delta_time,end_time, ...
        process_noise_variance_Q);
    measurement_array = measureState(state_array,measurement_noise_variance_R);
end

function state_x = generateState(init_state,start_time,delta_time, ...
        end_time,process_noise_variance_Q)
    
    time_eval_points = start_time:delta_time:end_time;
    state_x = zeros(1,numel(time_eval_points));
    state_x(1) = init_state;
    
    loop_run_counter = 0;
    for time = time_eval_points(2:end)
        loop_run_counter = loop_run_counter + 1;
        current_state = state_x(loop_run_counter);
        state_x(loop_run_counter+1) =  ...  % TODO: separate function
            current_state/2 + 25*current_state/(1 + current_state^2) ...
            + 8*cos(1.2*time) + normrnd(0,sqrt(process_noise_variance_Q)); 
    end
end

function measurement_z = measureState(state_x,measurement_noise_variance_R)
    state_vector_size = [1,numel(state_x)];
    measurement_z = state_x.^2./20 + ...
        normrnd(0,sqrt(measurement_noise_variance_R),state_vector_size);
end

function filtered_state = exactFlowFilter(particle_count,end_time,delta_time, ...
        process_noise_variance_Q,measurement_noise_variance_R,measurement_z)
    
    kLambdaDivPoints = 11;
    kInitialDistributionCovar = 200;
    
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
        state_array(:,1,i) = propagateState(state_array(:,1,i-1),process_noise_variance_Q,i);
        filtered_state(i) = mean(state_array(:,1,i));
        [predicted_EKF_mean_m(i),predicted_EKF_covar_P(i)] = ...
            EKFPrediction(EKF_mean_m(i-1),EKF_covar_P(i-1),process_noise_variance_Q,i);
        % iterate through the pseudo time lambda
        for j = 1:kLambdaDivPoints
            lambda = j * d_lambda;
            % iterate through the particles
            for k = 1:particle_count
                dx_dlambda = calculateSlope(lambda,state_array(k,j,i),filtered_state(i), ...
                    measurement_noise_variance_R,predicted_EKF_covar_P(i),measurement_z(i));
                state_array(k,j+1,i) = state_array(k,j,i) + dx_dlambda * d_lambda; 
            end
            filtered_state(i) = mean(state_array(:,j+1,i)); % re evaluate the avg
        end
        [EKF_mean_m(i),EKF_covar_P(i)] = EKFUpdate(predicted_EKF_mean_m(i), ...
            predicted_EKF_covar_P(i),measurement_noise_variance_R,measurement_z(i));
    end
end % TODO: get rid of the for loops -> vectorize

function dx_dlambda = calculateSlope(lambda,currentState,state_avg, ...
        measurement_noise_variance_R,predicted_EKF_covar_P,measurement_z)
    % Calculate the linearized measurement matrix for each particle
    % dx/dlambda = B*x + b
    
    linear_measurement_matrix_H = currentState./10; % TODO function
    B = -1/2*predicted_EKF_covar_P*linear_measurement_matrix_H' ...
        /(lambda*linear_measurement_matrix_H*predicted_EKF_covar_P ...
        *linear_measurement_matrix_H' + measurement_noise_variance_R) * linear_measurement_matrix_H;
    b = (eye(size(B)) + 2*lambda*B) ...
        * ((eye(size(B)) + lambda*B)*predicted_EKF_covar_P*linear_measurement_matrix_H'...
        / measurement_noise_variance_R*measurement_z + B*state_avg);
    
    dx_dlambda = B*currentState + b;
end

function propagated_state = propagateState(prev_state,process_noise_variance_Q,time)
    propagated_state = prev_state./2 + 25.*prev_state./(1 + prev_state.^2) ...
        + 8*cos(1.2*time) + normrnd(0,sqrt(process_noise_variance_Q),size(prev_state));
end

function [predicted_EKF_mean_m,predicted_EKF_covar_P] = EKFPrediction(EKF_mean_m,...
        EKF_covar_P,process_noise_variance_Q,time)
    
    predicted_EKF_mean_m = EKF_mean_m/2+25*EKF_mean_m/(1+EKF_mean_m^2) ...
        + 8*cos(1.2*time); % TODO function
    linear_system_matrix = 1/2 + 25 * (1 - EKF_covar_P^2) / ((1 + EKF_covar_P^2)^2); % Jacobian, TODO: function
    predicted_EKF_covar_P = linear_system_matrix*EKF_covar_P*linear_system_matrix' ...
        + process_noise_variance_Q;
end

function [EKF_mean_m,EKF_covar_P] = EKFUpdate(predicted_EKF_mean_m, ... 
        predicted_EKF_covar_P,measurement_noise_variance_R,measurement_z)
    
    linear_measurement_matrix_H = predicted_EKF_mean_m./10; % TODO function
    kalman_gain_K = predicted_EKF_covar_P * linear_measurement_matrix_H ...
        /(linear_measurement_matrix_H*predicted_EKF_covar_P*linear_measurement_matrix_H ...
        + measurement_noise_variance_R);
    EKF_mean_m = predicted_EKF_mean_m + kalman_gain_K ...
        *(measurement_z - expectedObservation(predicted_EKF_mean_m));
    EKF_covar_P = (1-kalman_gain_K*linear_measurement_matrix_H)*predicted_EKF_covar_P;
end

function z = expectedObservation(state)
    z = state.^2./20;  % TODO function
end
