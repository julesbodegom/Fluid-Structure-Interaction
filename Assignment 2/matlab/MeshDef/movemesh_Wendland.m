function nodesNew = movemesh_Wendland(nodes,disp,rSupport)
    
    Nb = sum(disp(:,1));       % boundary points (constraints)
    Ni = length(disp(:,1))-Nb; % internal points
    
    if (Ni == 0)
        % All points have prescribed displacements.
        nodesNew = nodes+disp(:,2:3);
        return;
    end

    Xb = zeros(Nb,2); % coordinates boundary points
    Db = zeros(Nb,2); % displacements boundary points
    Xi = zeros(Ni,2); % coordinates internal points
    Di = zeros(Ni,2); % displacements internal points
    ID = zeros(Ni,1); % global ID for internal points
    ib = 1;
    ii = 1;

    % Split nodes and disp into boundary and internal points.
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

    if (nargin < 3 || isempty(rSupport))
        rSupport = estimate_support_radius(Xb);
    end

    if (rSupport <= 0)
        error('movemesh_Wendland:InvalidSupportRadius', 'Support radius must be positive.');
    end

    % Build the compact-support Wendland C2 kernel block.
    Phi = zeros(Nb,Nb);
    for i = 1:Nb
        for j = 1:Nb
            rij = norm(Xb(i,:) - Xb(j,:));
            q = rij / rSupport;
            if (q <= 1)
                Phi(i,j) = (1-q)^4 * (4*q+1);
            else
                Phi(i,j) = 0;
            end
        end
    end

    % Linear polynomial augmentation (same structure as TPS implementation).
    P = [ones(Nb,1) Xb(:,1) Xb(:,2)];
    rbfMat = [Phi P;
              P' zeros(3,3)];

    rhs_x = [Db(:,1); zeros(3,1)];
    rhs_y = [Db(:,2); zeros(3,1)];

    coeff_x = rbfMat \ rhs_x;
    coeff_y = rbfMat \ rhs_y;

    gamma_x = coeff_x(1:Nb);
    poly_x  = coeff_x(Nb+1:Nb+3);
    gamma_y = coeff_y(1:Nb);
    poly_y  = coeff_y(Nb+1:Nb+3);

    % Evaluate displacement for each internal node.
    for i = 1:Ni
        x = Xi(i,1);
        y = Xi(i,2);

        dx_i = poly_x(1) + poly_x(2)*x + poly_x(3)*y;
        dy_i = poly_y(1) + poly_y(2)*x + poly_y(3)*y;

        for j = 1:Nb
            rij = norm(Xi(i,:) - Xb(j,:));
            q = rij / rSupport;
            if (q <= 1)
                phi = (1-q)^4 * (4*q+1);
            else
                phi = 0;
            end

            dx_i = dx_i + gamma_x(j) * phi;
            dy_i = dy_i + gamma_y(j) * phi;
        end

        Di(i,:) = [dx_i dy_i];
    end

    % Convert internal displacements to global disp vector.
    for i=1:Ni
        disp(ID(i),2:3) = Di(i,:);
    end

    nodesNew = nodes+disp(:,2:3);
end

function rSupport = estimate_support_radius(Xb)
    Nb = size(Xb,1);
    if (Nb < 2)
        rSupport = 1;
        return;
    end

    % Heuristic: three times the average nearest-neighbor spacing.
    nnDist = zeros(Nb,1);
    for i = 1:Nb
        dmin = inf;
        for j = 1:Nb
            if (i == j)
                continue;
            end
            dij = norm(Xb(i,:) - Xb(j,:));
            if (dij < dmin)
                dmin = dij;
            end
        end
        nnDist(i) = dmin;
    end

    rSupport = 3*mean(nnDist);

    if (~isfinite(rSupport) || rSupport <= 0)
        ext = max(Xb,[],1) - min(Xb,[],1);
        rSupport = 0.5*norm(ext);
    end
end