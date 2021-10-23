function[X_record, E_record] = MCMC_model_potential(kT, MCsteps, X_equilibrium, lower_boundary, upper_boundary, K, U)
%% Set the initial position to be "approximately" Target_X, e.g. -1.

U2 = @(x) U(x) + K.*(x - X_equilibrium).^2;
X_current  =  X_equilibrium;
E_current  =  U2(X_current); % compute the energy at the initial position.
sigma = 0.11; % this corresponds to the same deviation from currret position as the discrete metropolis with Num_X_points = 5 for 100 linearly spaced X_points.

X_record = zeros();
E_record = zeros();
%%
for i = 1:1:MCsteps
    %%
    X_new = X_current + sigma.*randn(1);
    
    if (X_new < lower_boundary || X_new > upper_boundary)
        
        X_record(i) = X_current;
        E_record(i) = E_current;
        
        continue
        
    end
    
    
    E_new = U2(X_new);
    dE = E_new - E_current;
    
    P = exp(-dE./kT);
    
    u = rand(1);
    
    if u <= min(1,P)
        
        E_current = E_new;
        X_current = X_new;
        
    end
    
    X_record(i) = X_current;
    E_record(i) = E_current;
end


end


