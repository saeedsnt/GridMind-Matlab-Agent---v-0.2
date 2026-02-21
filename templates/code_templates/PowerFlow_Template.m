%% Power Flow Analysis Template
% Auto-generated MATLAB code for load flow study
% 
% Modify parameters below and run to perform power flow analysis

clear; clc; close all;

%% System Parameters
baseMVA = 100;              % Base MVA
tolerance = 1e-6;           % Convergence tolerance
maxIter = 50;               % Maximum iterations
method = 'Newton-Raphson';  % Power flow method

%% Bus Data
% Format: [BusNum, Type, Vmag, Vang, Pg, Qg, Pl, Ql]
% Type: 1=PQ (Load), 2=PV (Generator), 3=Slack
busData = [
    1, 3, 1.000,  0.0,  0.00, 0.00, 0.00, 0.00;   % Slack bus
    2, 1, 1.000,  0.0,  0.00, 0.00, 50.0, 20.0;   % Load bus
    3, 2, 1.050,  0.0, 100.0, 0.00, 0.00, 0.00;   % Generator bus
    4, 1, 1.000,  0.0,  0.00, 0.00, 30.0, 10.0;   % Load bus
    5, 1, 1.000,  0.0,  0.00, 0.00, 40.0, 15.0;   % Load bus
];

%% Line Data
% Format: [FromBus, ToBus, R, X, B, Tap]
lineData = [
    1, 2, 0.0200, 0.0800, 0.0, 1.0;
    1, 3, 0.0150, 0.0600, 0.0, 1.0;
    2, 4, 0.0300, 0.1000, 0.0, 1.0;
    3, 4, 0.0250, 0.0850, 0.0, 1.0;
    4, 5, 0.0200, 0.0700, 0.0, 1.0;
];

%% Execute Power Flow
fprintf('\n================================================================================\n');
fprintf('             Power Flow Analysis - %s Method\n', method);
fprintf('================================================================================\n\n');

% Run Newton-Raphson or selected method
[V, converged] = runPowerFlow(busData, lineData, baseMVA, tolerance, maxIter);

if converged
    fprintf('Solution converged!\n\n');
    displayResults(busData, lineData, V, baseMVA);
else
    fprintf('Solution did not converge. Check system data.\n\n');
end

%% Function to run power flow
function [V, converged] = runPowerFlow(busData, lineData, baseMVA, tolerance, maxIter)
    n = size(busData,1);
    V = busData(:,3) .* exp(1i*busData(:,4)*pi/180);  % Initialize voltages
    
    % Build admittance matrix
    Y = zeros(n,n);
    for k = 1:size(lineData,1)
        f = lineData(k,1);
        t = lineData(k,2);
        R = lineData(k,3);
        X = lineData(k,4);
        B = lineData(k,5);
        tap = lineData(k,6);
        
        z = R + 1i*X;
        y = 1/z;
        
        Y(f,f) = Y(f,f) + y/tap^2 + 1i*B/2;
        Y(t,t) = Y(t,t) + y + 1i*B/2;
        Y(f,t) = Y(f,t) - y/tap;
        Y(t,f) = Y(t,f) - y/tap;
    end
    
    % Newton-Raphson iteration
    converged = false;
    for iter = 1:maxIter
        % Calculate power injections
        P = zeros(n,1);
        Q = zeros(n,1);
        for i = 1:n
            for j = 1:n
                P(i) = P(i) + real(V(i) * conj(V(j)) * Y(i,j));
                Q(i) = Q(i) + imag(V(i) * conj(V(j)) * Y(i,j));
            end
        end
        
        % Mismatch
        dP = (busData(:,5) - busData(:,7)) - P;
        dQ = (busData(:,6) - busData(:,8)) - Q;
        
        % Check convergence
        if max(abs([dP; dQ])) < tolerance
            converged = true;
            break;
        end
        
        % Update (simplified - add Jacobian for production code)
        V = V + 0.01 * (dP + 1i*dQ) ./ (abs(V) + 0.01);
    end
end

%% Display Results
function displayResults(busData, lineData, V, baseMVA)
    n = size(busData,1);
    
    fprintf('Bus   Voltage(pu)  Angle(deg)   Pgen(MW)  Pload(MW)  Qgen(MVAr) Qload(MVAr)\n');
    fprintf('-----  -----------  ---------   --------  ---------   ---------  ---------\n');
    
    for i = 1:n
        fprintf('%3d    %8.4f     %8.2f    %8.2f   %8.2f    %8.2f    %8.2f\n', ...
            i, abs(V(i)), angle(V(i))*180/pi, ...
            busData(i,5)*baseMVA, busData(i,7)*baseMVA, ...
            busData(i,6)*baseMVA, busData(i,8)*baseMVA);
    end
    
    % Line flows
    fprintf('\n\nLine Flows:\n');
    fprintf('From  To   Power(MW)  Reactive(MVAr)\n');
    fprintf('----  --   ---------  ---------------\n');
    
    for k = 1:size(lineData,1)
        f = lineData(k,1);
        t = lineData(k,2);
        R = lineData(k,3);
        X = lineData(k,4);
        
        S_f = V(f) * conj((V(f) - V(t))/(R + 1i*X));
        fprintf('%3d   %3d   %8.2f     %8.2f\n', f, t, real(S_f)*baseMVA, imag(S_f)*baseMVA);
    end
end
