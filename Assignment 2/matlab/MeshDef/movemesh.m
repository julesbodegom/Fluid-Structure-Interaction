function nodesNew = movemesh(nodes,disp)
    
    Nb = sum(disp(:,1));       % boundary points (constraints)
    Ni = length(disp(:,1))-Nb; % internal points
    
    if (Ni == 0)
        % All points have prescribed displacements
        % Update node locations
        nodesNew = nodes+disp(:,2:3);
    else
        Xb = zeros(Nb,2); % coordinates boundary points
        Db = zeros(Nb,2); % displacements boundary points
        Xi = zeros(Ni,2); % coordinates internal points
        Di = zeros(Ni,2); % displacements internal points
        ID = zeros(Ni,1); % global ID for internal points
        ib = 1;
        ii = 1;
        
        % Split nodes and disp into boundary and internal points
        for i=1:Nb+Ni
            if (disp(i,1) == 1)
                Xb(ib,:) = nodes(i,:);
                Db(ib,:) = disp(i,2:3);
                ib       = ib+1;
            else
                Xi(ii,:) = nodes(i,:);
                ID(ii)   = i;
                ii       = ii+1;
            end
        end
        
        % Set up RBF interpolation matrix for constraints
        rbfMat = zeros(Nb+3);

        %% IMPLEMENT YOUR RBF MESH DEFORMATION HERE
        % Build Phi block
        Phi = zeros(Nb,Nb);
        for i = 1:Nb
            for j = 1:Nb
                r = norm(Xb(i,:) - Xb(j,:));
                if r == 0
                    Phi(i,j) = 0;
                else
                    Phi(i,j) = r^2 * log(r);
                end
            end
        end

        % Build polynomial block
        P = [ones(Nb,1) Xb(:,1) Xb(:,2)];

        % Assemble full constraint matrix
        rbfMat = [Phi P;
                P' zeros(3,3)];
        

        
        % constraints = displacements of the boundary nodes (separate in x
        % and y)
        rhs_x = [Db(:,1); zeros(3,1)];
        rhs_y = [Db(:,2); zeros(3,1)];

        
        % determine rbf coefficients 
        coeff_x = rbfMat \ rhs_x;
        coeff_y = rbfMat \ rhs_y;

        
        % evaluate RBF function for each internal node
        gamma_x = coeff_x(1:Nb);
        poly_x  = coeff_x(Nb+1:Nb+3);

        gamma_y = coeff_y(1:Nb);
        poly_y  = coeff_y(Nb+1:Nb+3);
        % evaluate RBF function for each internal node
        for i = 1:Ni
            x = Xi(i,1);
            y = Xi(i,2);
            
            % polynomial part
            dx_i = poly_x(1) + poly_x(2)*x + poly_x(3)*y;
            dy_i = poly_y(1) + poly_y(2)*x + poly_y(3)*y;
            
            % RBF sum over all boundary points
            for j = 1:Nb
                r = norm(Xi(i,:) - Xb(j,:));
                if r ~= 0
                    phi = r^2 * log(r);
                else
                    phi = 0;
                end
                
                dx_i = dx_i + gamma_x(j) * phi;
                dy_i = dy_i + gamma_y(j) * phi;
            end
            
            Di(i,:) = [dx_i dy_i];
        end
                        
        % Convert internal displacements to global disp vector
        for i=1:Ni
            disp(ID(i),2:3) = Di(i,:);
        end
        
        % Update node locations
        nodesNew = nodes+disp(:,2:3);
    end
end