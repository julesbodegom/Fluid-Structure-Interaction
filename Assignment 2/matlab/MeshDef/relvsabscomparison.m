clear;
close all;
clc;

% Settings from the assignment
NdtP = 20;
N = 11;
M = 11;
maxAngle = 90;

% Run the four required cases
[t025_abs, q025_abs] = run_case(0.25, 1, NdtP, N, M, maxAngle);
[t025_rel, q025_rel] = run_case(0.25, 0, NdtP, N, M, maxAngle);

[t20_abs, q20_abs] = run_case(2.0, 1, NdtP, N, M, maxAngle);
[t20_rel, q20_rel] = run_case(2.0, 0, NdtP, N, M, maxAngle);

% Plot 1: Np = 0.25
figure(1);
clf;
plot(t025_abs, q025_abs, '-b', 'LineWidth', 1.5, 'DisplayName', 'Absolute displacement');
hold on;
plot(t025_rel, q025_rel, '-r', 'LineWidth', 1.5, 'DisplayName', 'Relative displacement');
grid on;
xlabel('time');
ylabel('min orthogonality');
title('Mesh Quality Comparison, Np = 0.25');
legend('Location', 'best');

% Plot 2: Np = 2.0
figure(2);
clf;
plot(t20_abs, q20_abs, '-g', 'LineWidth', 1.5, 'DisplayName', 'Absolute displacement');
hold on;
plot(t20_rel, q20_rel, '-c', 'LineWidth', 1.5, 'DisplayName', 'Relative displacement');
grid on;
xlabel('time');
ylabel('min orthogonality');
title('Mesh Quality Comparison, Np = 2.0');
legend('Location', 'best');

function [time, qualOfMesh] = run_case(Np, absdisp, NdtP, N, M, maxAngle)
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

    nSteps = round(NdtP * Np);
    qualOfMesh = ones(nSteps+1, 1);

    for iTimeStep = 1:nSteps
        t = iTimeStep * dt;
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
            nodes = movemesh(nodes0, displacements);
        else
            nodes = movemesh(nodes, displacements);
        end

        qualOfMesh(iTimeStep+1) = min_mesh_quality(cells, nodes);
    end

    time = (0:nSteps)' * dt;
end

function qmin = min_mesh_quality(mesh, nodes)
    qmin = 1;

    for iCell = 1:size(mesh, 1)
        Sn = compute_faceSn(mesh(iCell, :), nodes);
        n = nodes(mesh(iCell, :), :);
        q = 1;

        for iFace = 1:4
            v1 = Sn(iFace, :);
            v2 = n(mod(iFace+1,4)+1, :) - n(mod(iFace,4)+1, :);
            angle = (v1*v2') / sqrt((v1*v1') * (v2*v2'));
            q = min(q, (1-angle)/2);
        end

        qmin = min(qmin, q);
    end
end