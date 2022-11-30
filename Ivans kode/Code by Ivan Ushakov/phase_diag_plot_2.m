function phase_diag_plot_2(bg,bm,mats,coor)
% The function takes in spin liquid and gap limits, coordinate matrix
% (coor) stored in coordinates.mat and numerical matrix (mats), stored in
% solutions.mat. The function then plots the phase diagram. Remember to
% load the files where the data is stored, ad use the matrices as arguments
% in this function.
% bg: spin liquid limit (see eq. 5.4.81)
% bm: should be set to 0

[a,b] = size(mats);

col1 = [0.4940 0.1840 0.5560]; % sl
col2 = [0.3010 0.7450 0.9330]; % inc
col3 = [0.9290 0.6940 0.1250]; % com
col5 = [0 0 0];
col4 = [1 1 1];

% Begin plotting
ff = figure(01);
ss = 51;
plot(nan,nan,'s','MarkerSize',ss, 'MarkerEdgeColor',[0.8500 0.3250 0.0980],...
    'MarkerFaceColor',[0.8500 0.3250 0.0980])
hold on
drawnow;
plot(nan,nan,'s','MarkerSize',ss, 'MarkerEdgeColor',[0 0.4470 0.7410],...
    'MarkerFaceColor',[0 0.4470 0.7410])
plot(nan,nan,'s','MarkerSize',ss, 'MarkerEdgeColor',[0.4660 0.6740 0.1880],...
    'MarkerFaceColor',[0.4660 0.6740 0.1880])
plot(nan,nan,'s','MarkerSize',ss, 'MarkerEdgeColor',[0.4660 0.6740 0.1880],...
    'MarkerFaceColor',[1 1 1])
set(gca,'FontSize',20)
set(gcf,'Position',[100 100 600 600])
daspect([1,1,1])
set(gca,'FontSize',20)
xlabel('$\Gamma/J_1$','interpreter','Latex')
ylabel('$J_2/J_1$','interpreter','Latex')
% ------------------------------

for i = 1:a
    for j = 1:b
        
        % Phase diagram coordinates for the given step
        J2 = coor(3,i,j);
        g2 = coor(2,i,j);
        % ------------------------------
        
        % Solution matrix
        mat = mats{i,j};
        % ------------------------------
        
        vec = [];
        
        if ~isempty(mat)
            % Cut irrelevant solutions
            mat = mat_red(mat,4,0,1);       % all mu bigger than 0
            mat = mat_red(mat,6,1e-8,-1);   % small residue
            mat = mat_red(mat,9,-bm,1);      % gap restrictions
            % ------------------------------
            if ~isempty(mat)
                [~,aa] = min(mat(5,:));
                vec = mat(:,aa);
            end
        end
        
        if ~isempty(vec)
            if vec(9) > bg   % spin liquid
                plot(g2,J2,'s','MarkerSize',ss, 'MarkerEdgeColor',col1,...
    'MarkerFaceColor',col1)
                drawnow;    
            elseif abs(vec(8)) < 1e-2   % commensurate
                plot(g2,J2,'s','MarkerSize',ss, 'MarkerEdgeColor',col3,...
    'MarkerFaceColor',col3)
                drawnow;
            elseif vec(8) < 1.001*4*pi/(3*sqrt(3)) && ...
                    vec(8) > 0.999*pi/(sqrt(3))     % incommensurate
                plot(g2,J2,'s','MarkerSize',ss, 'MarkerEdgeColor',col2,...
    'MarkerFaceColor',col2)
                drawnow;
            else                                    % other
                plot(g2,J2,'s','MarkerSize',ss, 'MarkerEdgeColor',col5,...
    'MarkerFaceColor',col5)
                drawnow;
            end
        else                                    % no solution
            plot(g2,J2,'s','MarkerSize',ss, 'MarkerEdgeColor',col4,...
    'MarkerFaceColor',col4)
            drawnow;
        end
    end
end

gz = 0.4;
axis([-gz,6+gz,-gz,6+gz])
yticks(0:6)
xticks(0:6)
%set(gca,'visible','off')


end

