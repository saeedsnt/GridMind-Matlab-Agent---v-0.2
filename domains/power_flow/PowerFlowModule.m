classdef PowerFlowModule < handle
    % PowerFlowModule - Power flow analysis domain module
    %
    % Provides interactive parameter collection and code generation
    % for load flow studies using various methods.

    properties
        methodName          % Power flow method
        numBuses            % Number of buses
        tolerance           % Convergence tolerance
        maxIter             % Maximum iterations
        baseMVA             % Base MVA for per-unit system
        busData             % Bus data array
        lineData            % Line data array
    end

    methods
        function obj = PowerFlowModule()
            % Constructor - set defaults
            obj.methodName = 'Newton-Raphson';
            obj.numBuses = 5;
            obj.tolerance = 1e-6;
            obj.maxIter = 50;
            obj.baseMVA = 100;
            obj.busData = [];
            obj.lineData = [];
        end

        function params = collectParameters(obj, prompt)
            % Collect power flow parameters from user interactively

            fprintf('\n=== Power Flow Analysis Configuration ===\n\n');

            % Select method
            methods = {'Newton-Raphson', 'Gauss-Seidel', 'Fast Decoupled', 'DC Power Flow'};
            fprintf('Available power flow methods:\n');
            methodIdx = prompt.getMenuChoice(methods, 'Select power flow method:');
            obj.methodName = methods{methodIdx};
            fprintf('Selected method: %s\n\n', obj.methodName);

            % System size
            obj.numBuses = prompt.getNumericInputWithDefault( ...
                'Enter number of buses: ', 5, [2, 10000]);

            % Base MVA
            obj.baseMVA = prompt.getNumericInputWithDefault( ...
                'Enter base MVA: ', 100, [1, 10000]);

            % Convergence parameters
            obj.tolerance = prompt.getNumericInputWithDefault( ...
                'Enter convergence tolerance: ', 1e-6, [1e-12, 1e-3]);

            obj.maxIter = prompt.getNumericInputWithDefault( ...
                'Enter maximum iterations: ', 50, [10, 1000]);

            % Collect bus data
            obj.busData = obj.collectBusData(prompt, obj.numBuses);

            % Collect line data
            obj.lineData = obj.collectLineData(prompt, obj.numBuses);

            % Return parameters struct
            params.method = obj.methodName;
            params.numBuses = obj.numBuses;
            params.baseMVA = obj.baseMVA;
            params.tolerance = obj.tolerance;
            params.maxIter = obj.maxIter;
            params.busData = obj.busData;
            params.lineData = obj.lineData;
        end

        function busData = collectBusData(obj, prompt, numBuses)
            % Collect bus data from user
            busData = struct('type', {}, 'Vmag', {}, 'Vang', {}, 'Pg', {}, 'Qg', {}, 'Pl', {}, 'Ql', {});

            fprintf('\n--- Bus Data Entry ---\n');
            fprintf('Bus types: 1=PQ (Load), 2=PV (Generator), 3=Slack\n\n');

            for i = 1:numBuses
                fprintf('Bus %d:\n', i);

                % Bus type
                type = prompt.getNumericInput('  Bus type (1=PQ, 2=PV, 3=Slack): ', [1, 3]);
                busData(i).type = type;
                busData(i).name = sprintf('Bus%d', i);

                % Based on type, collect appropriate values
                switch type
                    case 1  % PQ bus
                        busData(i).Vmag = 1.0;
                        busData(i).Vang = 0;
                        busData(i).Pg = 0;
                        busData(i).Qg = 0;
                        busData(i).Pl = prompt.getNumericInputWithDefault('  Active load (MW): ', 50, [0, 10000]);
                        busData(i).Ql = prompt.getNumericInputWithDefault('  Reactive load (MVAr): ', 20, [0, 5000]);

                    case 2  % PV bus
                        busData(i).Vmag = prompt.getNumericInputWithDefault('  Voltage magnitude (pu): ', 1.05, [0.9, 1.1]);
                        busData(i).Vang = 0;
                        busData(i).Pg = prompt.getNumericInputWithDefault('  Active generation (MW): ', 100, [0, 10000]);
                        busData(i).Qg = 0;
                        busData(i).Pl = prompt.getNumericInputWithDefault('  Active load (MW): ', 0, [0, 10000]);
                        busData(i).Ql = prompt.getNumericInputWithDefault('  Reactive load (MVAr): ', 0, [0, 5000]);

                    case 3  % Slack bus
                        busData(i).Vmag = 1.0;
                        busData(i).Vang = 0;
                        busData(i).Pg = 0;  % To be determined
                        busData(i).Qg = 0;  % To be determined
                        busData(i).Pl = prompt.getNumericInputWithDefault('  Active load (MW): ', 0, [0, 10000]);
                        busData(i).Ql = prompt.getNumericInputWithDefault('  Reactive load (MVAr): ', 0, [0, 5000]);
                end
                fprintf('\n');
            end
        end

        function lineData = collectLineData(obj, prompt, numBuses)
            % Collect line/branch data from user
            lineData = struct('fromBus', {}, 'toBus', {}, 'R', {}, 'X', {}, 'B', {}, 'tap', {});

            fprintf('\n--- Line Data Entry ---\n');
            fprintf('Enter transmission line/transformer data\n\n');

            numLines = prompt.getNumericInput('Enter number of lines/transformers: ', [1, numBuses*(numBuses-1)]);

            for i = 1:numLines
                fprintf('Line %d:\n', i);

                % From and To buses
                fromBus = prompt.getNumericInput('  From bus number: ', [1, numBuses]);
                toBus = prompt.getNumericInput('  To bus number: ', [1, numBuses]);

                while fromBus == toBus
                    fprintf('  From and To buses must be different.\n');
                    toBus = prompt.getNumericInput('  To bus number: ', [1, numBuses]);
                end

                lineData(i).fromBus = fromBus;
                lineData(i).toBus = toBus;
                lineData(i).name = sprintf('Line%d', i);

                % Line parameters (in per-unit)
                lineData(i).R = prompt.getNumericInputWithDefault('  Resistance R (pu): ', 0.02, [0, 1]);
                lineData(i).X = prompt.getNumericInputWithDefault('  Reactance X (pu): ', 0.08, [0, 1]);
                lineData(i).B = prompt.getNumericInputWithDefault('  Shunt susceptance B (pu): ', 0, [0, 1]);

                % Transformer tap (1.0 for transmission lines)
                isTransformer = prompt.getYesNo('  Is this a transformer? (y/n): ');
                if isTransformer
                    lineData(i).tap = prompt.getNumericInputWithDefault('  Tap ratio: ', 1.0, [0.8, 1.2]);
                else
                    lineData(i).tap = 1.0;
                end
                fprintf('\n');
            end
        end

        function code = generateCode(obj, params)
            % Generate MATLAB code for power flow analysis
            code = obj.generatePowerFlowScript();
        end

        function script = generatePowerFlowScript(obj)
            % Generate complete power flow script

            script = [
                '%% Power Flow Analysis' newline
                '% Auto-generated MATLAB code for load flow study' newline
                newline
                'clear; clc; close all;' newline
                newline
                '%% System Parameters' newline
                sprintf('baseMVA = %g;' newline, obj.baseMVA) newline
                sprintf('tolerance = %g;' newline, obj.tolerance) newline
                sprintf('maxIter = %d;' newline, obj.maxIter) newline
                sprintf('method = '''', obj.methodName)
                ''';' newline
                newline
                '%% Bus Data' newline
                '% Format: [BusNum, Type, Vmag, Vang, Pg, Qg, Pl, Ql]' newline
                '% Type: 1=PQ, 2=PV, 3=Slack' newline
                'busData = [' newline
            ];

            % Add bus data
            for i = 1:length(obj.busData)
                bus = obj.busData(i);
                script = [script, sprintf('    %d, %d, %.3f, %.1f, %.2f, %.2f, %.2f, %.2f', ...
                    i, bus.type, bus.Vmag, bus.Vang, bus.Pg/baseMVA, bus.Qg/baseMVA, ...
                    bus.Pl/baseMVA, bus.Ql/baseMVA)];
                if i < length(obj.busData)
                    script = [script, ';' newline];
                else
                    script = [script, ';' newline];
                end
            end
            script = [script, '];' newline newline];

            % Add line data
            script = [script
                '%% Line Data' newline
                '% Format: [FromBus, ToBus, R, X, B, Tap]' newline
                'lineData = [' newline
            ];

            for i = 1:length(obj.lineData)
                line = obj.lineData(i);
                script = [script, sprintf('    %d, %d, %.4f, %.4f, %.4f, %.3f', ...
                    line.fromBus, line.toBus, line.R, line.X, line.B, line.tap)];
                if i < length(obj.lineData)
                    script = [script, ';' newline];
                else
                    script = [script, ';' newline];
                end
            end
            script = [script, '];' newline newline];

            % Add power flow function based on method
            switch obj.methodName
                case 'Newton-Raphson'
                    script = [script, obj.getNewtonRaphsonCode()];
                case 'Gauss-Seidel'
                    script = [script, obj.getGaussSeidelCode()];
                case 'Fast Decoupled'
                    script = [script, obj.getFastDecoupledCode()];
                case 'DC Power Flow'
                    script = [script, obj.getDCPowerFlowCode()];
            end

            % Add results display
            script = [script
                newline
                '%% Display Results' newline
                'fprintf(''\n========================================\n'')' newline
                'fprintf(''       POWER FLOW RESULTS\n'')' newline
                'fprintf(''========================================\n\n'')' newline
                'fprintf(''Bus   V(pu)    Angle(deg)   Pgen    Pload    Qgen    Qload\n'')' newline
                'fprintf(''------------------------------------------------------------\n'')' newline
                'for i = 1:size(busData,1)' newline
                '    fprintf(''%3d   %7.4f   %9.3f     %6.3f   %6.3f   %6.3f   %6.3f\n'', ...' newline
                '        i, V(i)*180/pi*1i, angle(V(i))*180/pi, ...' newline
                '        real(Sg(i)), real(Sl(i)), imag(Sg(i)), imag(Sl(i)));' newline
                'end' newline
                'fprintf(''\nTotal Generation: %.3f + j%.3f pu\n'', sum(real(Sg)), sum(imag(Sg)));' newline
                'fprintf(''Total Load:       %.3f + j%.3f pu\n'', sum(real(Sl)), sum(imag(Sl)));' newline
                'fprintf(''Total Losses:     %.3f + j%.3f pu\n'', sum(real(Sg))-sum(real(Sl)), ...' newline
                '    sum(imag(Sg))-sum(imag(Sl)));' newline
                'fprintf(''========================================\n'')' newline
            ];

            % Add optional plotting
            script = [script
                newline
                '%% Voltage Profile Plot' newline
                'figure;' newline
                'bar(1:length(V), abs(V));' newline
                'xlabel(''Bus Number'');' newline
                'ylabel(''Voltage Magnitude (pu)'');' newline
                'title(''Voltage Profile'');' newline
                'grid on;' newline
                'ylim([0.9 1.1]);' newline
                'yline(1.0, ''k--'');' newline
            ];
        end

        function code = getNewtonRaphsonCode(obj)
            % Generate Newton-Raphson power flow code
            code = [
                '%% Newton-Raphson Power Flow' newline
                'n = size(busData,1);' newline
                'V = ones(n,1);  % Initialize voltages' newline
                'theta = zeros(n,1);' newline
                '' newline
                '% Build admittance matrix' newline
                'Y = zeros(n,n);' newline
                'for k = 1:size(lineData,1)' newline
                '    f = lineData(k,1);' newline
                '    t = lineData(k,2);' newline
                '    R = lineData(k,3);' newline
                '    X = lineData(k,4);' newline
                '    B = lineData(k,5);' newline
                '    tap = lineData(k,6);' newline
                '    z = R + 1i*X;' newline
                '    y = 1/z;' newline
                '    Y(f,f) = Y(f,f) + y/tap^2 + 1i*B/2;' newline
                '    Y(t,t) = Y(t,t) + y + 1i*B/2;' newline
                '    Y(f,t) = Y(f,t) - y/tap;' newline
                '    Y(t,f) = Y(t,f) - y/tap;' newline
                'end' newline
                '' newline
                '% Net power injection (generation - load)' newline
                'Pspec = (busData(:,5) - busData(:,7)) * baseMVA;' newline
                'Qspec = (busData(:,6) - busData(:,8)) * baseMVA;' newline
                'busType = busData(:,2);' newline
                '' newline
                '% Newton-Raphson iteration' newline
                'for iter = 1:maxIter' newline
                '    % Calculate power injections' newline
                '    Pcalc = zeros(n,1);' newline
                '    Qcalc = zeros(n,1);' newline
                '    for i = 1:n' newline
                '        for j = 1:n' newline
                '            Pcalc(i) = Pcalc(i) + abs(V(i))*abs(V(j))*... ' newline
                '                abs(Y(i,j))*cos(angle(Y(i,j)) - theta(i) + theta(j));' newline
                '            Qcalc(i) = Qcalc(i) - abs(V(i))*abs(V(j))*... ' newline
                '                abs(Y(i,j))*sin(angle(Y(i,j)) - theta(i) + theta(j));' newline
                '        end' newline
                '    end' newline
                '    ' newline
                '    % Mismatch calculation' newline
                '    dP = Pspec - Pcalc;' newline
                '    dQ = Qspec - Qcalc;' newline
                '    ' newline
                '    % Check convergence for PQ and PV buses' newline
                '    pqBuses = find(busType == 1);' newline
                '    pvBuses = find(busType == 2);' newline
                '    slackBuses = find(busType == 3);' newline
                '    ' newline
                '    mismatchP = dP([pqBuses; pvBuses]);' newline
                '    mismatchQ = dQ(pqBuses);' newline
                '    mismatch = [mismatchP; mismatchQ];' newline
                '    ' newline
                '    if max(abs(mismatch)) < tolerance' newline
                '        fprintf(''Converged in %d iterations\n'', iter);' newline
                '        break;' newline
                '    end' newline
                '    ' newline
                '    % Build Jacobian matrix (simplified - full implementation needed)' newline
                '    % For production use, implement full Jacobian calculation' newline
                '    ' newline
                '    % Update angles and voltages (simplified)' newline
                '    theta([pqBuses; pvBuses]) = theta([pqBuses; pvBuses]) + 0.1*dP([pqBuses; pvBuses]);' newline
                '    V(pqBuses) = V(pqBuses) + 0.1*dQ(pqBuses);' newline
                'end' newline
                '' newline
                '% Convert to complex voltages' newline
                'V = V .* exp(1i*theta);' newline
                '' newline
                '% Calculate final power flows' newline
                'Sg = (busData(:,5) + 1i*busData(:,6));' newline
                'Sl = (busData(:,7) + 1i*busData(:,8));' newline
            ];
        end

        function code = getGaussSeidelCode(obj)
            % Generate Gauss-Seidel power flow code
            code = [
                '%% Gauss-Seidel Power Flow' newline
                'n = size(busData,1);' newline
                'V = ones(n,1);  % Flat start' newline
                'accelFactor = 1.6;  % Acceleration factor' newline
                '' newline
                '% Build admittance matrix' newline
                'Y = zeros(n,n);' newline
                'for k = 1:size(lineData,1)' newline
                '    f = lineData(k,1);' newline
                '    t = lineData(k,2);' newline
                '    R = lineData(k,3);' newline
                '    X = lineData(k,4);' newline
                '    B = lineData(k,5);' newline
                '    tap = lineData(k,6);' newline
                '    z = R + 1i*X;' newline
                '    y = 1/z;' newline
                '    Y(f,f) = Y(f,f) + y/tap^2 + 1i*B/2;' newline
                '    Y(t,t) = Y(t,t) + y + 1i*B/2;' newline
                '    Y(f,t) = Y(f,t) - y/tap;' newline
                '    Y(t,f) = Y(t,f) - y/tap;' newline
                'end' newline
                '' newline
                '% Net power (generation - load) in per-unit' newline
                'Pnet = (busData(:,5) - busData(:,7));' newline
                'Qnet = (busData(:,6) - busData(:,8));' newline
                'busType = busData(:,2);' newline
                '' newline
                '% Gauss-Seidel iteration' newline
                'for iter = 1:maxIter' newline
                '    Vold = V;' newline
                '    ' newline
                '    for i = 1:n' newline
                '        if busType(i) == 3  % Slack bus' newline
                '            continue;' newline
                '        end' newline
                '        ' newline
                '        % Calculate new voltage' newline
                '        sumYV = 0;' newline
                '        for j = 1:n' newline
                '            if j ~= i' newline
                '                sumYV = sumYV + Y(i,j) * V(j);' newline
                '            end' newline
                '        end' newline
                '        ' newline
                '        Snet = Pnet(i) - 1i*Qnet(i);' newline
                '        Vnew = (conj(Snet)/conj(V(i)) - sumYV) / Y(i,i);' newline
                '        ' newline
                '        % Apply acceleration' newline
                '        V(i) = V(i) + accelFactor * (Vnew - V(i));' newline
                '        ' newline
                '        % For PV buses, maintain voltage magnitude' newline
                '        if busType(i) == 2' newline
                '            V(i) = busData(i,3) * V(i) / abs(V(i));' newline
                '        end' newline
                '    end' newline
                '    ' newline
                '    % Check convergence' newline
                '    if max(abs(V - Vold)) < tolerance' newline
                '        fprintf(''Converged in %d iterations\n'', iter);' newline
                '        break;' newline
                '    end' newline
                'end' newline
                '' newline
                '% Calculate power flows' newline
                'Sg = (busData(:,5) + 1i*busData(:,6));' newline
                'Sl = (busData(:,7) + 1i*busData(:,8));' newline
            ];
        end

        function code = getFastDecoupledCode(obj)
            % Generate Fast Decoupled power flow code
            code = [
                '%% Fast Decoupled Power Flow (XB Version)' newline
                'n = size(busData,1);' newline
                'V = ones(n,1);' newline
                'theta = zeros(n,1);' newline
                '' newline
                '% Build B'' and B'' matrices (simplified)' newline
                '% For production use, implement full fast decoupled formulation' newline
                'Bprime = zeros(n,n);  % For P-theta' newline
                'Bdoubleprime = zeros(n,n);  % For Q-V' newline
                '' newline
                'for k = 1:size(lineData,1)' newline
                '    f = lineData(k,1);' newline
                '    t = lineData(k,2);' newline
                '    X = lineData(k,4);' newline
                '    b = -1/X;  % Susceptance' newline
                '    Bprime(f,t) = b;' newline
                '    Bprime(t,f) = b;' newline
                '    Bprime(f,f) = Bprime(f,f) - b;' newline
                '    Bprime(t,t) = Bprime(t,t) - b;' newline
                '    Bdoubleprime(f,t) = b;' newline
                '    Bdoubleprime(t,f) = b;' newline
                '    Bdoubleprime(f,f) = Bdoubleprime(f,f) - b;' newline
                '    Bdoubleprime(t,t) = Bdoubleprime(t,t) - b;' newline
                'end' newline
                '' newline
                'Pnet = (busData(:,5) - busData(:,7));' newline
                'Qnet = (busData(:,6) - busData(:,8));' newline
                'busType = busData(:,2);' newline
                '' newline
                'for iter = 1:maxIter' newline
                '    % Calculate power mismatch' newline
                '    dP = Pnet;' newline
                '    dQ = Qnet;' newline
                '    ' newline
                '    % Solve decoupled equations (simplified)' newline
                '    pqBuses = find(busType == 1);' newline
                '    pvBuses = find(busType == 2);' newline
                '    ' newline
                '    % Update angles' newline
                '    theta(pqBuses) = theta(pqBuses) + 0.1*dP(pqBuses);' newline
                '    theta(pvBuses) = theta(pvBuses) + 0.1*dP(pvBuses);' newline
                '    ' newline
                '    % Update voltages' newline
                '    V(pqBuses) = V(pqBuses) + 0.1*dQ(pqBuses);' newline
                '    ' newline
                '    if max(abs([dP; dQ])) < tolerance' newline
                '        fprintf(''Converged in %d iterations\n'', iter);' newline
                '        break;' newline
                '    end' newline
                'end' newline
                '' newline
                'V = V .* exp(1i*theta);' newline
                'Sg = (busData(:,5) + 1i*busData(:,6));' newline
                'Sl = (busData(:,7) + 1i*busData(:,8));' newline
            ];
        end

        function code = getDCPowerFlowCode(obj)
            % Generate DC power flow code
            code = [
                '%% DC Power Flow' newline
                'n = size(busData,1);' newline
                '' newline
                '% Build B matrix (susceptance only)' newline
                'B = zeros(n,n);' newline
                'for k = 1:size(lineData,1)' newline
                '    f = lineData(k,1);' newline
                '    t = lineData(k,2);' newline
                '    X = lineData(k,4);' newline
                '    b = -1/X;' newline
                '    B(f,t) = b;' newline
                '    B(t,f) = b;' newline
                '    B(f,f) = B(f,f) - b;' newline
                '    B(t,t) = B(t,t) - b;' newline
                'end' newline
                '' newline
                '% Net power injection (MW)' newline
                'P = (busData(:,5) - busData(:,7));' newline
                '' newline
                '% Remove slack bus row/column (assume bus 1 is slack)' newline
                'Bred = B(2:end, 2:end);' newline
                'Pred = P(2:end);' newline
                '' newline
                '% Solve for angles' newline
                'theta_red = Bred \ Pred;' newline
                'theta = [0; theta_red];  % Slack bus angle = 0' newline
                '' newline
                '% Calculate line flows' newline
                'fprintf(''DC Power Flow Results:\n'');' newline
                'for k = 1:size(lineData,1)' newline
                '    f = lineData(k,1);' newline
                '    t = lineData(k,2);' newline
                '    X = lineData(k,4);' newline
                '    Pflow = (theta(f) - theta(t)) / X;' newline
                '    fprintf(''Line %d->%d: %.3f pu (%.1f MW)\n'', f, t, Pflow, Pflow*baseMVA);' newline
                'end' newline
                '' newline
                'V = ones(n,1) .* exp(1i*theta);  % DC approximation' newline
                'Sg = (busData(:,5) + 1i*busData(:,6));' newline
                'Sl = (busData(:,7) + 1i*busData(:,8));' newline
            ];
        end

        function generateSimulinkModel(obj, params)
            % Generate Simulink model for power flow visualization
            fprintf('\nSimulink model generation for power flow:\n');
            fprintf('Note: Power flow is typically a code-based calculation.\n');
            fprintf('For dynamic studies, consider using the stability module.\n\n');

            % Generate model creation script
            modelScript = obj.generateSimulinkScript();
            fprintf('%s\n', modelScript);
        end

        function script = generateSimulinkScript(obj)
            % Generate script to create Simulink model
            script = [
                '% Power Flow Visualization Model' newline
                'modelName = ''power_flow_viz'';' newline
                'if system_exists(modelName)' newline
                '    close_system(modelName, 0);' newline
                '    delete_system(modelName);' newline
                'end' newline
                'new_system(modelName);' newline
                '' newline
                '% Add blocks for bus visualization' newline
                '% (Implementation requires Simscape Electrical toolbox)' newline
                'fprintf(''Model created: %s.slx\n'', modelName);' newline
            ];
        end

        function handleInput(obj, userInput)
            % Handle contextual user input
            fprintf('Power flow module received: %s\n', userInput);
        end
    end
end
