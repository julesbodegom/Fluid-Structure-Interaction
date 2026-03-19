clear;
close all;
clc;

% Settings from previous run
NdtP = 20;
Np   = 0.25;
N = 11;
M = 11;
maxAngle = 90;

% We'll use absolute displacement to clearly show the deformation
absdisp = 0; 

radii = [0.25, 0.75];

for iRad = 1:length(radii)
    rSupport = radii(iRad);
    
    [cells, nodes] = get_final_mesh(Np, absdisp, NdtP, N, M, maxAngle, rSupport);
    
    figure(iRad);
    clf;
    plotMesh(cells, nodes);
    colorbar;
    title(sprintf('Final Mesh Quality (r = %.2f, relative displacement)', rSupport));
    xlabel('x');
    ylabel('y');
end

function [cells, nodes] = get_final_mesh(Np, absdisp, NdtP, N, M, maxAngle, rSupport)
    Lx = 1;
    Ly = 1;

    % Initial mesh
    nodes = zeros((N+1)*(M+1), 2);
    for x = 0:M
        nodes(x*(N+1)+1:(x+1)*(N+1), 1) = x/M*Lx;
        nodes(x*(N+1)+1:(x+1)*(N+1), 2) = (0:N)'/N*Ly;
    end

    cells = create_mesh(N, M, nodes);

    omega = 2*pi;
    dt = 1/NdtP;
    rc = [Lx/2, Ly/2];
    nodes0 = nodes;

    % Moving nodes
    id1 = (N+1)*floor((M-1)/2) + floor((N+1)/2);
    id2 = id1 + 1;
    id3 = id1 + (N+1);
    id4 = id3 + 1;
    dispIDs = [id1 id2 id3 id4];

    % Static boundary nodes
    staticIDs = [1:N+1, (1:M-1)*(N+1)+1, (2:M)*(N+1), M*(N+1)+1:(M+1)*(N+1)];

    nSteps = round(NdtP*Np);

    for iTimeStep = 1:nSteps
        t = iTimeStep*dt;
        theta_old = pi/180 * maxAngle * sin(omega*(t-dt));
        theta = pi/180 * maxAngle * sin(omega*t);

        displacements = zeros((N+1)*(M+1), 3);

        for nodeID = dispIDs
            if absdisp == 1
                xold = nodes0(nodeID, 1);
                yold = nodes0(nodeID, 2);
                dTheta = theta;
            else
                xold = nodes(nodeID, 1);
                yold = nodes(nodeID, 2);
                dTheta = theta - theta_old;
            end

            xnew = rc(1) + (xold-rc(1))*cos(dTheta) - (yold-rc(2))*sin(dTheta);
            ynew = rc(2) + (yold-rc(2))*cos(dTheta) + (xold-rc(1))*sin(dTheta);

            displacements(nodeID, :) = [1, xnew-xold, ynew-yold];
        end

        for nodeID = staticIDs
            displacements(nodeID, :) = [1, 0, 0];
        end

        if absdisp == 1
            nodes = movemesh_Wendland(nodes0, displacements, rSupport);
        else
            nodes = movemesh_Wendland(nodes, displacements, rSupport);
        end
    end
end
