function di = deltak(kx, ky, D1, D2, D3, J1, J2, G, D)

%   deltak
%   The function finds the term inside the squareroot of E^{-}_k in 
%   eq. 5.4.78, except for mu^2. This is used to compute the gap.

       gk = J1^2*D1^2/4*(1 + 4*cos(sqrt(3)/2*kx)*cos(3/2*ky) + ...          % This value corresponds to |\tilde{\eta}_k|^2 in eq. 5.4.56, and is rewritten in terms of real functions only for the sake of numerical stability
           4*(cos(sqrt(3)/2*kx))^2);
       
       Sk = 4*sin(sqrt(3)/2*kx)*(cos(3/2*ky) - cos(sqrt(3)/2*kx));          % eta_k eq. 5.4.58
       
       Ck = 2*(2*cos(sqrt(3)/2*kx)*cos(3/2*ky) + cos(sqrt(3)*kx));          % zeta_k eq. 5.4.59
       
       es = 0.5*J2*D2*Sk;                                                   % Im(psi_k) in eq. 5.4.56
       
       et = 0.5*G*D3*Ck;                                                    % Re(psi_k) in eq. 5.4.56
       
       zs = 0.25*D*D2*Ck;                                                   % Im(tau_k) in eq. 5.4.56
       
       zt = 0.25*D*D3*Sk;                                                   % Re(tau_k) in eq. 5.4.56
       
       gz = sqrt(gk*(zt.^2 + es.^2) + (et.*zt + es.*zs).^2);            % This value corresponds to the squareroot of eq. 5.4.79 divided by 2.
       
       
       di = - gk - es.^2 - et.^2 - zs.^2 - zt.^2 - 2*gz;                 % This value corresponds to all the terms except mu^2 inside the squareroot of E^{-}_k in eq. 5.4.78. Minimizing this function brings us to condensation point
       

end

