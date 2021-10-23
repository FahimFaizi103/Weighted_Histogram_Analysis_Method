
%%
clear all
clc

kT = 1.2;
MCsteps = 10^4;

%% Model potential.

G = 10;
U = @(x) G.*(x -1).^2.*(x + 1).^2;

figure
hold on
fplot(U,[-1.5, 1.5],'k','Linewidth',2)

%%
lower_boundary =  -1.5;
upper_boundary = 1.5;
sigma = 0.11;

num_windows = 14;
X_equilibrium = linspace(lower_boundary, upper_boundary, num_windows);

X_equilibrium_record = zeros();

K = ones(1,num_windows).*200;
%%
X_record = zeros(num_windows, MCsteps);
E_record = zeros(num_windows, MCsteps);

for i = 1:1:num_windows
    
   
    [X_record(i,:), E_record(i,:)] = MCMC_model_potential(kT, MCsteps, X_equilibrium(i), lower_boundary, upper_boundary, K(i), U);    
    X_equilibrium_record(i) = X_equilibrium(i);
    
    
end

save('umbrella_sampling_data.mat')
