function F = sumzz_bos(df, x, J1, J2, G, D, N)

%   sumzz_bos 
%   selfconsistency equations to be solved
%   This function takes the desired self-consistency functions as argument,
%   the initial guess, and the numerical parameters. Returns k-sums of the
%   self-consistency equations, which for a solution should be equal to 0

F = zeros(1, length(df));

for i = 1:length(df)
   F(i) = ksumm_bos(df{i}, x(1), x(2), x(3), x(4), J1, J2, G, D, N); 
end

end