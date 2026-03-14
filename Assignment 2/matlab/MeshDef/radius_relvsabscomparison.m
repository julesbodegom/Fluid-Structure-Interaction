clear;
close all;
clc;

% Requested settings
NdtP = 20;
N = 11;
M = 11;
maxAngle = 90;

% Simulation horizon (same style as relvsabscomparison)
Np = 2.0;

% Four Wendland support radii to investigate
radii = [0.2 0.4 0.6 0.8 1.0];
nr = length(radii);

% Store quality histories for each radius (columns).
qAbsAll = [];
qRelAll = [];
t = [];

for iRad = 1:nr
    rSupport = radii(iRad);

    [t_abs, q_abs] = run_case(Np, 1, NdtP, N, M, maxAngle, rSupport);
    [t_rel, q_rel] = run_case(Np, 0, NdtP, N, M, maxAngle, rSupport);

    if isempty(t)
        t = t_abs;
        qAbsAll = zeros(length(q_abs), nr);
        qRelAll = zeros(length(q_rel), nr);
    end

    qAbsAll(:, iRad) = q_abs;
    qRelAll(:, iRad) = q_rel;
end

% Plot 1: absolute displacement, all radii.
figure(1);
clf;
hold on;
for iRad = 1:nr
    plot(t, qAbsAll(:, iRad), 'LineWidth', 1.5, 'DisplayName', sprintf('r = %.2f', radii(iRad)));
end
grid on;
xlabel('time');
ylabel('min orthogonality');
title('Absolute displacement: support-radius comparison');
legend('Location', 'best');

% Plot 2: relative displacement, all radii.
figure(2);
clf;
hold on;
for iRad = 1:nr
    plot(t, qRelAll(:, iRad), 'LineWidth', 1.5, 'DisplayName', sprintf('r = %.2f', radii(iRad)));
end
grid on;
xlabel('time');
ylabel('min orthogonality');
title('Relative displacement: support-radius comparison');
legend('Location', 'best');

% Additional plot: final min orthogonality vs support radius
% Requested settings: NdtP = 20, Np = 0.25
NpRadius = 0.25;
radiiSweep = [0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0];
nSweep = length(radiiSweep);
finalAbs = zeros(nSweep, 1);
finalRel = zeros(nSweep, 1);

for iRad = 1:nSweep
    rSupport = radiiSweep(iRad);

    [~, q_abs] = run_case(NpRadius, 1, NdtP, N, M, maxAngle, rSupport);
    [~, q_rel] = run_case(NpRadius, 0, NdtP, N, M, maxAngle, rSupport);

    finalAbs(iRad) = q_abs(end);
    finalRel(iRad) = q_rel(end);
end

figure(3);
clf;
plot(radiiSweep, finalAbs, '-ob', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Absolute displacement');
hold on;
plot(radiiSweep, finalRel, '-sr', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Relative displacement');
grid on;
xlabel('support radius');
ylabel('final min orthogonality');
title('Final min orthogonality vs support radius (Np = 0.25, NdtP = 20)');
legend('Location', 'best');

function [time, qualOfMesh] = run_case(Np, absdisp, NdtP, N, M, maxAngle, rSupport)
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
    qualOfMesh = ones(nSteps+1, 1);

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