function[bias] = compute_bias(K, x, X_equilibrium)

bias = K.*(x - X_equilibrium).^2;


end