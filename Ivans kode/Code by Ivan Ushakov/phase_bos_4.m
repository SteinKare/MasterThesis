function M = phase_bos_4(df, J1, J2, G, D, N, x0)

%   phase_bos_4 
%   The function returns essential information about soultions of
%   self-consistency equations
%   The function takes selfconsistency k-terms, interaction parameters,
%   lattice size and initial guesses as arguments. The function solves the
%   four derivatives of eq. 5.4.76 = 0. The function returns a matrix
%   that contains the mean-field parameters that satisfy the equations,
%   success number and numerical residue (see below)


% solve self-consistency equations with fsolve
fun = @ (x) sumzz_bos(df, x, J1, J2, G, D, N);                              % self-consistency equations with k-summation
options = optimset('Display','off');                                        % turn off display
[x,fval,exitflag] = fsolve(fun,x0,options);                                 % the solution
% --------------------------

% construct the solution matrix
M = zeros(1,6);

M(1, 1) = x(1);                                                             % computed D1                          
M(1, 2) = x(2);                                                             % computed D2                     
M(1, 3) = x(3);                                                             % computed D3            
M(1, 4) = x(4);                                                             % computed mu
M(1, 5) = exitflag;                                                         % solution success number
M(1, 6) = abs(fval(1)) + abs(fval(2)) + ...
    abs(fval(3)) + abs(fval(4));                                            % numerical residue
% --------------------------


end

