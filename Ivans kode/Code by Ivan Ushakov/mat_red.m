function mat_red = mat_red(mat,par,lim,p)
%   mat_red
%   Used to eliminate irrelevant solutions
%   The function takes in a matrix (mat), matrix row of interest (par), 
%   limit value (lim) and a value p (1 or -1). The function then removes 
%   all the matrix columns where the row value number par is smaller
%   (p = 1) or greater (p = -1) than lim value.

if p == 1
    logic = real(mat(par,:)) > lim;
elseif p == -1
    logic = real(mat(par,:)) < lim;
end

S = sum(logic);
[l,~] = size(mat);
mat_red = zeros(l,S);

for j = 1:l
   A = mat(j,:);
   mat_red(j,:) = A(logic);
end

end


