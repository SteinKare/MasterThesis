function S = energy_bos_ref(kx, ky, D1, D2, D3, mu, J1, J2, G, D, kappa)

% the function takes in the k-vector, mean-field parameters and interaction
% parameters. The function returns the single term in eq. 5.4.76 for D = 0.
% The function is only used as a test and reference for D = 0.
       
       gk = J1^2*D1^2/4*(1 + 4*cos(sqrt(3)/2*kx)*cos(3/2*ky) + ...          % this value corresponds to |\tilde{\eta}_k|^2 in eq. 5.4.56, and is rewritten in terms of real functions only for the sake of numerical stability
           4*(cos(sqrt(3)/2*kx))^2);
       
       es = 2*J2*D2*sin(sqrt(3)/2*kx)*(cos(3/2*ky) - cos(sqrt(3)/2*kx));    % Im(\psi_k) in eq. 5.4.56 
       
       et = G*D3*(2*cos(sqrt(3)/2*kx)*cos(3/2*ky) + cos(sqrt(3)*kx));       % Re(\psi_k) in eq. 5.4.56
       
       Ep = sqrt(mu^2 + 2*sqrt(gk)*abs(es) - gk - ...                       % eq. 5.4.77 (plus version)
           es^2 - et^2);
       
       Em = sqrt(mu^2 - 2*sqrt(gk)*abs(es) - gk - ...                       % eq. 5.4.77 (minus version)
           es^2 - et^2);
       
       S = Ep + Em + 3/2*J1*D1^2 + 3*J2*D2^2 + 3*G*D3^2 - (2+2*kappa)*mu;   % the sum term in eq. 5.4.76

end

