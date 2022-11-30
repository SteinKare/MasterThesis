function S = ksum(fun, N)
%	ksum 
%   First Brillouin zone k-sum 
%   The function takes in a function (of k) as argument, and then returns
%   the sum (per site) of this function in the first Brilluoin zone (see
%   proposition 5.5)
%   In this script, kx = kx*a and ky = ky*a

S = 0;
for nx = 1:N
    for ny = 1:N
       
       kx = (nx - ny)*2*pi/(N*sqrt(3));
       ky = (nx + ny)*2*pi/(3*N);
       kx = kx + 1e-6;                                                      % without this, the calculations became unstable
       ky = ky + 1e-6;                                                      % without this, the calculations became unstable
       M = fun(kx, ky);
       
       S = S + M;
    end
    
end

S = S/(N^2); % gives sum per site.

end

