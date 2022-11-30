function [val,kx] = ktopp(D1,D2,D3,J1,J2,G,D)
%   ktop 
%   Function that determines dispersion minima
%   The function takes in the three mean-field parameters and interaction 
%   parameters, and returns the greatest value (in k-space in the two 
%   disperision) of the quantity in the squareroot of the dispersion that 
%   is substracted from \mu (x_4). The function also returns the given 
%   k-value


ky = 0;                                                                     % we search for minima along kx_axis, as it contains both Gamma- and K-point (where we expect to find minima)
z = 4*pi/(sqrt(3)*3);                                                       % K-point

su = @(kx) deltak(kx, ky, D1, D2, D3, J1, J2, G, D);                        
x01 = z*(1 + 3/4)/2;                                                        % checking minima near K-point
x02 = 0;                                                                    % checking minima near Gamma-point
options = optimset('Display','off');
[x1,fval1,exitflag1] = fminsearch(su,x01,options);
[x2,fval2,exitflag2] = fminsearch(su,x02,options);

if exitflag1 ~= 1
    fval1 = input('No solution found. Enter fval1 manually: ');
    x1 = input('No solution found. Enter x1 manually: ');
end

if exitflag2 ~= 1
    fval2 = input('No solution found. Enter fval2 manually: ');
    x2 = input('No solution found. Enter x2 manually: ');
end

if fval1 < fval2
    val = -fval1;
    kx = x1;
else
    val = -fval2;
    kx = x2;
end

