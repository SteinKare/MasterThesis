function [df,f] = symbolic_diff(kappa,settings)

% The function takes in boson density (kappa) and settings parameter (see
% explanations in energy_bos functions). The function returns the k-terms
% in eq. 5.4.76 (stored in f) and it's partial derivatives with respect to D1, D2, D3
% and mu (stored in df).

    % symbolic variables (used when differentiating)
    syms kx_s ky_s D1_s D2_s D3_s mu_s J1_s J2_s G_s D_s
    % ---------------------
    
    % energy functions (see explanations inside the functions)
    
    if isequal(settings, 'no absolute value')
        f(kx_s, ky_s, D1_s, D2_s, D3_s, mu_s, J1_s, J2_s, G_s, D_s) = ...
            energy_bos_ref(kx_s, ky_s, D1_s, D2_s, ...
            D3_s, mu_s, J1_s, J2_s, G_s, D_s, kappa);
        
    elseif isequal(settings, 'absolute value')
        f(kx_s, ky_s, D1_s, D2_s, D3_s, mu_s, J1_s, J2_s, G_s, D_s) = ...
            energy_bos(kx_s, ky_s, D1_s, D2_s, ...
            D3_s, mu_s, J1_s, J2_s, G_s, D_s, kappa);
        
    end
    f = simplify(f);                                                        % f (sum term in eq. 5.4.76) is symbolically simplified
    % ---------------------

    % first self-consistency equation
    dif1 = diff(f, D1_s);                                                   % f is differentiated with respect to Delta_1
    dif1 = simplify(dif1);                                                  % differentiated function is symbolically simplified
    dif1 = matlabFunction(dif1);                                            % differentiated function is converted from symbolical to matlab function
    % ---------------------
    
    % second self-consistency equation
    dif2 = diff(f, D2_s);                                                   % f is differentiated with respect to Delta_2
    dif2 = simplify(dif2);
    dif2 = matlabFunction(dif2);
    % ---------------------
    
    % third self-consistency equation
    dif3 = diff(f, D3_s);                                                   % f is differentiated with respect to Delta_3
    dif3 = simplify(dif3);
    dif3 = matlabFunction(dif3);
    % ---------------------
    
    % forth self-consistency equation
    dif4 = diff(f, mu_s);                                                   % f is differentiated with respect to mu
    dif4 = simplify(dif4);
    dif4 = matlabFunction(dif4);
    % ---------------------

    df = {dif1, dif2, dif3, dif4};                                          % differentiated functions are stored in a structure array
    f = matlabFunction(f);                                                  % non-differentiated function

end

