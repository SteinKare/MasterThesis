function S = energy_bos(kx, ky, D1, D2, D3, mu, J1, J2, G, D, kappa)

% the function takes in the k-vector, mean-field parameters and interaction
% parameters. The function returns the single term in eq. 5.4.76 for any D.

       gk = J1^2*D1^2/4*(1 + 4*cos(sqrt(3)/2*kx)*cos(3/2*ky) + ...          % This value corresponds to |\tilde{\eta}_k|^2 in eq. 5.4.56, and is rewritten in terms of real functions only for the sake of numerical stability
           4*(cos(sqrt(3)/2*kx))^2);
       
       Sk = 4*sin(sqrt(3)/2*kx)*(cos(3/2*ky) - cos(sqrt(3)/2*kx));          % eta_k eq. 5.4.58
       
       Ck = 2*(2*cos(sqrt(3)/2*kx)*cos(3/2*ky) + cos(sqrt(3)*kx));          % zeta_k eq. 5.4.59
       
       es = 0.5*J2*D2*Sk;                                                   % Im(psi_k) in eq. 5.4.56
       
       et = 0.5*G*D3*Ck;                                                    % Re(psi_k) in eq. 5.4.56
       
       zs = 0.25*D*D2*Ck;                                                   % Im(tau_k) in eq. 5.4.56
       
       zt = 0.25*D*D3*Sk;                                                   % Re(tau_k) in eq. 5.4.56
       
       gz = sqrt(gk*(zt.^2 + es.^2) + (et.*zt + es.*zs).^2);            % This value corresponds to the squareroot of eq. 5.4.79 divided by 2.
       
       Ep = sqrt(abs(mu^2 - gk - es.^2 - et.^2 - zs.^2 - zt.^2 + 2*gz)); % plus version of eq. 5.4.80
       
       Em = sqrt(abs(mu^2 - gk - es.^2 - et.^2 - zs.^2 - zt.^2 - 2*gz)); % minus version of eq. 5.4.80
       
       % eq. 5.4.80 without absolute value transformation:
       
       %Ep = sqrt(mu^2 + 2*sqrt(gk)*abs(es) - gk - ...
       %    es^2 - et^2);
       %Em = sqrt(mu^2 - 2*sqrt(gk)*abs(es) - gk - ...
       %    es^2 - et^2);
       
       % ---------------------------------------------
       
       S = Ep + Em + 3/2*J1*D1^2 + 3*J2*D2^2 + 3*G*D3^2 - (2+2*kappa)*mu;

end

