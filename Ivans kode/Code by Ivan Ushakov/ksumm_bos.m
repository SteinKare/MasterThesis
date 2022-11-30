function S = ksumm_bos(func, D1, D2, D3, mu, J1, J2, G, D, N)

%   ksumm 
%   k-sum of a more specific function
%   This function combines a desired function that is to be summed in
%   k-space with relevant numerical values

fun = @ (kx, ky) func(kx, ky, D1, D2, D3, mu, J1, J2, G, D);
S = ksum(fun, N);

end

