function S = phase_solver(J1,J2,G,D,N,u,df,f)
%   phase_solver 
%   Solution of self-consistency equations
%   The function takes in physical paramters, as well as initial guesses,
%   total energy k-terms and self-consistency k-terms

% defining solution matrices
uuu = [];                                                                   % here we store indices of initial guesses where solution is found
mat = [];                                                                   % here we store the solutions and other relevant information (see below)
% --------------------------

% finding number of initial guesses
[~,b] = size(u);
% --------------------------


for i = 1:b
    % progress information
    %clc
    %fprintf('Point %d out of %d\n', i,b)
    % --------------------------
    
    % using phase_bos_4 to find solutions and ksumm_bos to find total energy
    x0 = [u(1,i),u(2,i),u(3,i),u(4,i)];                                     % initial guess
    x = phase_bos_4(df, J1, J2*J1, G*J1, D, N, x0);                         % self-consistency eq. solutions
    x = real(x);                                                            % sometimes the solutions enter imaginary domain
    
    %{
    % experiment
    k = 1;
    eeee = @ (kx, ky, D1s, D2s, D2t, u, J1, J2, g2, D) energy_bos_4(kx, ky, D1s, D2s, D2t, u, J1, J2, g2, D, k);
    e = ksumm_bos(eeee, x(1), x(2), x(3), x(4), J1, J2, g2, D, N);
    e = real(e);
    % --------------------------
    %}
    
    e = ksumm_bos(f, x(1), x(2), x(3), x(4), J1, J2, G, D, N);              % Finding total energy for our parameters (eq. 5.4.76)
    f1 = ksumm_bos(df{1}, x(1), x(2), x(3), x(4), J1, J2, G, D, N);         % Self-consistency eq. 1 (should be approx. 0)
    f2 = ksumm_bos(df{2}, x(1), x(2), x(3), x(4), J1, J2, G, D, N);         % Self-consistency eq. 2 (should be approx. 0)
    f3 = ksumm_bos(df{3}, x(1), x(2), x(3), x(4), J1, J2, G, D, N);         % Self-consistency eq. 3 (should be approx. 0)
    f4 = ksumm_bos(df{4}, x(1), x(2), x(3), x(4), J1, J2, G, D, N);         % Self-consistency eq. 4 (should be approx. 0)
    % --------------------------
    
    % putting results into solution matrices
    if x(5) >= 1                                                            % Only if solution is found
        uuu(1,end+1) = i;                                                   % initial guess index
        mat(1,end+1) = x(1);                                                % Computed D1
        mat(2,end) = x(2);                                                  % Computed D2
        mat(3,end) = x(3);                                                  % Computed D3
        mat(4,end) = x(4);                                                  % Computed mu
        mat(5,end) = e;                                                     % Resulting total energy
        mat(6,end) = abs(f1) + abs(f2) + abs(f3) + abs(f4);                 % residue
        
        % finding relevant k-point and gap size
        [valq,kxq] = ktopp(x(1),x(2),x(3),J1,J2,G,D);
        mat(7,end) = x(4)^2 - valq;                                         % Inside the square root value
        disp("Inside the sqrt value:")
        disp(x(4)^2 - valq)
        mat(8,end) = kxq;                                                   % Condensation/gap point
        % --------------------------
    end
    % --------------------------
    
end

if ~isempty(mat)
    mat(9,:) = mat(7,:)./(mat(4,:).^2);                                     % relative gap (see eq. 5.4.81)
    disp("Gap")
        disp(mat(9,:))
end

S = {uuu,mat};

end



