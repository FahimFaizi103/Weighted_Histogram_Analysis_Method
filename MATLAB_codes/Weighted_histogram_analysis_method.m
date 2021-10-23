%% References

% 1. Extensions to the weighted histogram analysis method: combining umbrella sampling with free energy calculations, Souaille and Roux
% 2. Catalytic Mechanism of RNA Backbone Cleavage by Ribonuclease H from QM/MM Simulations by Rosta et.al

%% Load the data

clear all
clc

load('umbrella_sampling_data.mat')

%% Compute Biased Probability density

Num_Bins = 100;
window = linspace(1,num_windows,num_windows);

figure
hold on
for i = 1:1:num_windows
    
    [values,edges]=histcounts(X_record(window(i),:),Num_Bins);
    values = values/sum(values);
    
    plot(edges(2:end), values, 'Linewidth', 2) % probability distribution
%     plot(edges(2:end), -log(values)-min(-log(values)), 'Linewidth',2) % free energy
    
end

%% Computing F_i values

beta = 1./kT;
diff_threshold = 10^(-3);
diff = 10*diff_threshold;

F_record = [];
F = ones(1,num_windows);
F_record(1,:) = F;
%%
tt = 1;
while diff > diff_threshold
    
    tt = tt+1;
    
    bias = zeros(num_windows, MCsteps);
    numerator = zeros(num_windows, MCsteps, num_windows);
    denominator = zeros(num_windows, MCsteps);
    
    for kk = 1:1:num_windows
        
        bias = compute_bias(K(kk), X_record, X_equilibrium(kk));
        numerator(:,:,kk) = exp(-beta.*bias);
        
        denominator = denominator + MCsteps.*exp(-beta.*bias + beta.*F(kk));
        
    end
    
    arg = numerator./denominator;
    arg = sum(sum(arg,2));
    
    F = -kT.*log(arg(:));
    F_record = [F_record; F'];
    
    diff = abs(max(F_record(tt,:) - F_record(tt-1,:)));
    
    
end

%% plot the F_i values against num_iterations to visualize convergence

 figure
 plot(F_record)
 
 %% Computing the unbiased probability distribution.
 
Num_Bins = 100;
EDGES = linspace(lower_boundary,upper_boundary,Num_Bins+1);
beta = 1./kT;
prob = zeros();
for i = 2:1:length(EDGES)

    if i == length(EDGES)
        
    candidates = find(X_record >= EDGES(i-1) & X_record <= EDGES(i) );
    
    else
        
    candidates = find(X_record >= EDGES(i-1) & X_record < EDGES(i) );
        
    end
    
    Sum = 0;
    for kk = 1:1:length(candidates)
        
        denom = 0;
        for j = 1:1:num_windows
            
            bias = compute_bias(K(j), EDGES(i), X_equilibrium(j));
            denom = denom + MCsteps.*exp(-beta.*bias + beta.*F(j));
                      
        end
        
        Sum = Sum + 1./denom;
              
    end
    
    prob(i-1) = Sum;
    
end

prob = prob./sum(prob);

%% Plot the probability distribution
figure
plot(EDGES(2:end), prob, 'Linewidth', 2) % probability distribution

%% Plot the free energy profile
figure
hold on
plot(EDGES(2:end), -log(prob)-min(-log(prob)), 'Linewidth',2)

% True Free energy plot

hold on;
G = 10;
beta = 1./kT;
U = @(x) beta.*G.*(x -1).^2.*(x + 1).^2; % times by beta, because in your code below you are computing free energy profile PER Temperature. i.e. F/kT = -log(values)
fplot(U,[-1.5, 1.5],'k', 'Linewidth',2)

legend('WHAM','Exact')
