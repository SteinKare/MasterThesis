tic

% Doing symbolig differentiation
kappa = 0.3139;                                                             % Boson density (see eq. 5.4.19)
settings = 'absolute value';                                                % See explanation in symbolic_diff
[df,f] = symbolic_diff(kappa,settings);                                     % f is the k-terms in eq. 5.4.76; df are its derivatives w.r.t. mean-field parameters
% ------------------------------

% Begin plotting
figure(01)                                                                  % The phase diagram is plotted as it is calculated
plot(nan, nan, 'r*', 'linewidth', 5)
hold on
drawnow;
plot(nan, nan, 'b*', 'linewidth', 5)
plot(nan, nan, 'g*', 'linewidth', 5)
plot(nan, nan, 'c*', 'linewidth', 5)
set(gca,'FontSize',20)
daspect([1,1,1])
xlabel('$\Gamma/J_1$','interpreter','Latex')
ylabel('$J_2/J_1$','interpreter','Latex')
col1 = [0.4940 0.1840 0.5560]; % spin liquid color
col2 = [0.3010 0.7450 0.9330]; % incommensurte order color
col3 = [0.9290 0.6940 0.1250]; % commensurate order color
col5 = [0 0 0];                % other solution color
col4 = [0.9 0.9 0.9];          % no solution color
% ------------------------------

% Physical and computation parameters
J1 = 1;                                                                     % NN interaction is set to 1, as all other interactions can be divided by J1
D = 0;                                                                 % DMI interaction
N = 20;                                                                     % number of lattice sites in each sublattice in one spatial direction
l = 200;                                                                     % number of random initial guesses for each phase diagram point
b = 0.005;                                                                  % spin liquid limit (see eq. 5.4.81)  

u = rand(4,l);                                                              % Random initial mean-field parameters (D1, D2, D3, mu), each between 0-1
u = ones(4,l);    
    % Scaling of random initial parameters 
    u = u*1;                                                                % Scaling D1, D2, D3 with 1
    u(4,:) = u(4,:)*10;                                                     % Scaling mu with 10
    % ------------------------------
% ------------------------------


% Phase diagram parameters
Z = 10;                                                                     % Number of phase diagram points in each direction
dlim = 6;                                                                   % phase diagram limit for J2 and Gamma
Gv = linspace(0,dlim,Z);                                                    % Vector of G values
J2v = linspace(0,dlim,Z);                                                   % Vector of J2 values
Gv(1) = 1e-9;                                                               % Values equal to exactly 0 tend to give unstable solutions, so we set very small values instead
J2v(1) = 1e-9;
% ------------------------------

% More plot settings
ss = 275/Z;                                                                 
gz = dlim/15;
axis([-gz,dlim+gz,-gz,dlim+gz])
yticks(0:dlim)
xticks(0:dlim)
% ------------------------------

% Load initial file (phase diagram data is saved here)
inp = input('Do you want to delete old phase diagram data? (1 = yes, 0 = no): ');
load('coordinates.mat')
load('solutions.mat')
if inp == 1
    coor = [];
    mats = [];
end
% ------------------------------

for i = 1:length(Gv)
    for j = 1:length(J2v)
        
        % Phase diagram coordinates for the given step
        J2 = J2v(j);
        G = Gv(i);
        % ------------------------------
        
        % Save coordinates to matrix
        coor(2,i,j) = G;
        coor(3,i,j) = J2;
        % ------------------------------
        
        % Solution matrix
        S = phase_solver(J1,J2,G,D,N,u,df,f);
        mat = S{2};
        mats{i,j} = mat;
        % ------------------------------
        
        vec = [];
        
        if ~isempty(mat)
            
            % Cut irrelevant solutions
            mat = mat_red(mat,4,0,1);                                       % all mu bigger than 0
            mat = mat_red(mat,6,1e-8,-1);                                   % small residue
            mat = mat_red(mat,9,0,1);                                       % avoid imaginary solutions (which is the case if mat(9, -) is less than 0)
            % ------------------------------
            
            % Pick the solution with lowest total energy
            if ~isempty(mat)                                                
                [~,a] = min(mat(5,:));
                vec = mat(:,a);
            end
            % ------------------------------

        end
        
        % Plot phase diagram
        if ~isempty(vec)
            if vec(9) > b   % spin liquid
                plot(G,J2,'s','MarkerSize',ss, 'MarkerEdgeColor',col1,...
                'MarkerFaceColor',col1)
                drawnow;    
                coor(1,i,j) = 2;
            elseif abs(vec(8)) < 1e-2   % commensurate (kx near Gamma)
                plot(G,J2,'s','MarkerSize',ss, 'MarkerEdgeColor',col3,...
                'MarkerFaceColor',col3)
                drawnow;
                coor(1,i,j) = 1;
            elseif vec(8) < 1.001*4*pi/(3*sqrt(3)) && ...
                    vec(8) > 0.999*pi/(sqrt(3))     % incommensurate (kx near K)
                plot(G,J2,'s','MarkerSize',ss, 'MarkerEdgeColor',col2,...
                'MarkerFaceColor',col2)
                drawnow;
                coor(1,i,j) = 0;
            else                                    % other
                plot(G,J2,'s','MarkerSize',ss, 'MarkerEdgeColor',col5,...
                'MarkerFaceColor',col5)
                drawnow;
                coor(1,i,j) = -2;
            end
        else                                    % no solution
            plot(G,J2,'s','MarkerSize',ss, 'MarkerEdgeColor',col4,...
            'MarkerFaceColor',col4)
            drawnow;
            coor(1,i,j) = -1;
        end
        % ------------------------------
        
        % Save the data to a file
        save('coordinates','coor')      % coordinates file
        save('solutions','mats')      % numerical results matrix
        % ------------------------------
        
    end
end

% finalizing plot
gz = 10*dlim/(15*Z);
axis([-gz,dlim+gz,-gz,dlim+gz])
yticks(0:dlim)
xticks(0:dlim)
% -------------------------

t = toc;