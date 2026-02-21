classdef StabilityModule < handle
    % StabilityModule - Stability analysis domain module
    %
    % Provides interactive parameter collection and code generation
    % for power system stability studies including transient,
    % small signal, and voltage stability analysis.

    properties
        studyType             % Type of stability study
        simulationTime        % Simulation duration (seconds)
        stepSize              % Time step for simulation
        solver                % ODE solver type
        generatorData         % Generator parameters
        disturbanceType       % Type of disturbance
        faultClearingTime     % Fault clearing time
        loadData              % Load parameters
        networkData           % Network parameters
    end

    methods
        function obj = StabilityModule()
            % Constructor - set defaults
            obj.studyType = 'transient';
            obj.simulationTime = 10;
            obj.stepSize = 0.001;
            obj.solver = 'ode23t';
            obj.generatorData = [];
            obj.disturbanceType = 'three_phase_fault';
            obj.faultClearingTime = 0.1;
            obj.loadData = [];
            obj.networkData = [];
        end

        function params = collectParameters(obj, prompt)
            % Collect stability analysis parameters from user

            fprintf('\n=== Stability Analysis Configuration ===\n\n');

            % Select study type
            studyTypes = {'transient', 'small_signal', 'voltage', 'frequency'};
            studyDescriptions = {
                'Transient Stability (time-domain simulation)', ...
                'Small Signal Stability (eigenvalue analysis)', ...
                'Voltage Stability (P-V/Q-V curves)', ...
                'Frequency Stability (load-frequency dynamics)'
            };

            fprintf('Available study types:\n');
            studyIdx = prompt.getSelectOne(studyTypes, studyDescriptions, 'Select study type:');
            obj.studyType = studyTypes{studyIdx};
            fprintf('Selected study type: %s\n\n', obj.studyType);

            % Simulation parameters
            obj.simulationTime = prompt.getNumericInputWithDefault( ...
                'Enter simulation time (seconds): ', 10, [1, 100]);

            obj.stepSize = prompt.getNumericInputWithDefault( ...
                'Enter time step (seconds): ', 0.001, [1e-5, 0.1]);

            % Solver selection
            solvers = {'ode45', 'ode23t', 'ode15s', 'ode23'};
            solverIdx = prompt.getMenuChoice(solvers, 'Select ODE solver:');
            obj.solver = solvers{solverIdx};

            % Disturbance type
            disturbances = {
                'three_phase_fault', 'line_trip', 'generator_trip', ...
                'load_step', 'voltage_sag'
            };
            distIdx = prompt.getMenuChoice(disturbances, 'Select disturbance type:');
            obj.disturbanceType = disturbances{distIdx};

            if contains(obj.disturbanceType, 'fault')
                obj.faultClearingTime = prompt.getNumericInputWithDefault( ...
                    'Enter fault clearing time (seconds): ', 0.1, [0.01, 1]);
            end

            % Collect system data
            obj.generatorData = obj.collectGeneratorData(prompt);
            obj.loadData = obj.collectLoadData(prompt);
            obj.networkData = obj.collectNetworkData(prompt);

            % Return parameters struct
            params.studyType = obj.studyType;
            params.simulationTime = obj.simulationTime;
            params.stepSize = obj.stepSize;
            params.solver = obj.solver;
            params.disturbanceType = obj.disturbanceType;
            params.faultClearingTime = obj.faultClearingTime;
            params.generatorData = obj.generatorData;
            params.loadData = obj.loadData;
            params.networkData = obj.networkData;
        end

        function genData = collectGeneratorData(obj, prompt)
            % Collect generator data for stability study
            genData = struct('bus', {}, 'M', {}, 'D', {}, 'Xd', {}, 'Xdp', {}, ...
                             'Tdo', {'Ra', 'H', 'type'});

            fprintf('\n--- Generator Data ---\n');
            numGens = prompt.getNumericInput('Enter number of generators: ', [1, 50]);

            for i = 1:numGens
                fprintf('Generator %d:\n', i);
                genData(i).bus = prompt.getNumericInput('  Connected bus: ', [1, 100]);
                genData(i).H = prompt.getNumericInputWithDefault( ...
                    '  Inertia constant H (seconds): ', 4.0, [1, 10]);
                genData(i).D = prompt.getNumericInputWithDefault( ...
                    '  Damping coefficient D (pu): ', 0.01, [0, 1]);
                genData(i).Xd = prompt.getNumericInputWithDefault( ...
                    '  Synchronous reactance Xd (pu): ', 1.8, [0.5, 5]);
                genData(i).Xdp = prompt.getNumericInputWithDefault( ...
                    '  Transient reactance X''d (pu): ', 0.3, [0.1, 1]);
                genData(i).Tdo = prompt.getNumericInputWithDefault( ...
                    '  Open circuit time constant T''d0 (seconds): ', 8.0, [1, 20]);
                genData(i).Ra = prompt.getNumericInputWithDefault( ...
                    '  Armature resistance Ra (pu): ', 0.003, [0, 0.1]);
                genData(i).Pg = prompt.getNumericInputWithDefault( ...
                    '  Initial active power (pu): ', 0.8, [0, 2]);
                genData(i).Vt = prompt.getNumericInputWithDefault( ...
                    '  Terminal voltage (pu): ', 1.05, [0.9, 1.1]);
                fprintf('\n');
            end
        end

        function loadData = collectLoadData(obj, prompt)
            % Collect load data for stability study
            loadData = struct('bus', {}, 'P', {}, 'Q', {}, 'type', {});

            fprintf('\n--- Load Data ---\n');
            numLoads = prompt.getNumericInput('Enter number of loads: ', [1, 100]);

            for i = 1:numLoads
                fprintf('Load %d:\n', i);
                loadData(i).bus = prompt.getNumericInput('  Connected bus: ', [1, 100]);
                loadData(i).P = prompt.getNumericInputWithDefault( ...
                    '  Active power (pu): ', 0.5, [0, 10]);
                loadData(i).Q = prompt.getNumericInputWithDefault( ...
                    '  Reactive power (pu): ', 0.2, [0, 5]);

                % Load type
                loadTypes = {'constant_impedance', 'constant_power', 'constant_current', 'composite'};
                typeIdx = prompt.getMenuChoice(loadTypes, '  Load model type:');
                loadData(i).type = loadTypes{typeIdx};
                fprintf('\n');
            end
        end

        function netData = collectNetworkData(obj, prompt)
            % Collect network data for stability study
            netData = struct('fromBus', {}, 'toBus', {}, 'R', {}, 'X', {}, 'B', {});

            fprintf('\n--- Network Data ---\n');
            numLines = prompt.getNumericInput('Enter number of lines: ', [1, 200]);

            for i = 1:numLines
                fprintf('Line %d:\n', i);
                netData(i).fromBus = prompt.getNumericInput('  From bus: ', [1, 100]);
                netData(i).toBus = prompt.getNumericInput('  To bus: ', [1, 100]);
                netData(i).R = prompt.getNumericInputWithDefault( ...
                    '  Resistance R (pu): ', 0.01, [0, 1]);
                netData(i).X = prompt.getNumericInputWithDefault( ...
                    '  Reactance X (pu): ', 0.05, [0, 1]);
                netData(i).B = prompt.getNumericInputWithDefault( ...
                    '  Shunt susceptance B (pu): ', 0, [0, 1]);
                fprintf('\n');
            end
        end

        function code = generateCode(obj, params)
            % Generate MATLAB code for stability analysis
            switch obj.studyType
                case 'transient'
                    code = obj.generateTransientStabilityCode();
                case 'small_signal'
                    code = obj.generateSmallSignalCode();
                case 'voltage'
                    code = obj.generateVoltageStabilityCode();
                case 'frequency'
                    code = obj.generateFrequencyStabilityCode();
                otherwise
                    code = obj.generateTransientStabilityCode();
            end
        end

        function script = generateTransientStabilityCode(obj)
            % Generate transient stability analysis code

            script = [
                '%% Transient Stability Analysis' newline
                '% Auto-generated MATLAB code for transient stability study' newline
                newline
                'clear; clc; close all;' newline
                newline
                '%% System Parameters' newline
                sprintf('f = 60;  %% System frequency (Hz)' newline)
                sprintf('omega_s = 2*pi*f;  %% Synchronous speed (rad/s)' newline)
                sprintf('simTime = %g;  %% Simulation time (s)' newline, obj.simulationTime)
                sprintf('dt = %g;  %% Time step (s)' newline, obj.stepSize)
                sprintf('tFault = 0.1;  %% Fault application time (s)' newline)
                sprintf('tClear = %g;  %% Fault clearing time (s)' newline, obj.faultClearingTime)
                newline
                '%% Generator Parameters (Classical Model)' newline
            ];

            % Add generator parameters
            if isempty(obj.generatorData)
                script = [script
                    '% Single machine infinite bus (SMIB) example' newline
                    'H = 4.0;  %% Inertia constant (s)' newline
                    'D = 0.01;  %% Damping coefficient (pu)' newline
                    'Xd_prime = 0.3;  %% Transient reactance (pu)' newline
                    'Pm = 0.8;  %% Mechanical power input (pu)' newline
                    'Vt = 1.05;  %% Terminal voltage (pu)' newline
                    'Vinf = 1.0;  %% Infinite bus voltage (pu)' newline
                    'Xline = 0.5;  %% Line reactance (pu)' newline
                ];
            else
                gen = obj.generatorData(1);
                script = [script
                    sprintf('H = %g;  %% Inertia constant (s)\n', gen.H)
                    sprintf('D = %g;  %% Damping coefficient (pu)\n', gen.D)
                    sprintf('Xd_prime = %g;  %% Transient reactance (pu)\n', gen.Xdp)
                    sprintf('Pm = %g;  %% Mechanical power input (pu)\n', gen.Pg)
                    sprintf('Vt = %g;  %% Terminal voltage (pu)\n', gen.Vt)
                    'Vinf = 1.0;  %% Infinite bus voltage (pu)' newline
                    'Xline = 0.5;  %% Line reactance (pu)' newline
                ];
            end

            script = [script
                newline
                '%% Initial Conditions' newline
                '% Calculate initial rotor angle from power flow' newline
                'Xtotal = Xd_prime + Xline;' newline
                'delta0 = asin(Pm * Xtotal / (Vt * Vinf));  %% Initial rotor angle' newline
                'omega0 = 0;  %% Initial speed deviation' newline
                newline
                '%% Swing Equation Parameters' newline
                'M = 2 * H / omega_s;  %% Inertia coefficient' newline
                'Pmax = Vt * Vinf / Xtotal;  %% Maximum power transfer' newline
                newline
                '%% Time Vector' newline
                't = 0:dt:simTime;' newline
                'nSteps = length(t);' newline
                newline
                '%% Pre-allocate Arrays' newline
                'delta = zeros(1, nSteps);' newline
                'omega = zeros(1, nSteps);' newline
                'delta(1) = delta0;' newline
                'omega(1) = omega0;' newline
                newline
                '%% Fault Modeling' newline
                '% During fault: Pmax drops significantly' newline
                'Pmax_prefault = Pmax;' newline
                'Pmax_during_fault = 0.1 * Pmax;  %% Fault reduces transfer capability' newline
                'Pmax_post_fault = Pmax;  %% Post-fault (line may be out)' newline
                newline
                '%% Simulation Loop (Euler Method)' newline
                'for i = 2:nSteps' newline
                '    % Determine system condition' newline
                '    if t(i) < tFault' newline
                '        Pmax_current = Pmax_prefault;' newline
                '    elseif t(i) < tClear' newline
                '        Pmax_current = Pmax_during_fault;' newline
                '    else' newline
                '        Pmax_current = Pmax_post_fault;' newline
                '    end' newline
                '    ' newline
                '    % Swing equations' newline
                '    Pe = Pmax_current * sin(delta(i-1));  %% Electrical power output' newline
                '    domega_dt = (Pm - Pe - D * omega(i-1)) / M;  %% Acceleration' newline
                '    ddelta_dt = omega(i-1);  %% Rate of change of angle' newline
                '    ' newline
                '    % Euler integration' newline
                '    omega(i) = omega(i-1) + domega_dt * dt;' newline
                '    delta(i) = delta(i-1) + ddelta_dt * dt;' newline
                'end' newline
                newline
            ];

            % Add results and plots
            script = [script
                '%% Display Results' newline
                'fprintf(''\n========================================\n'')' newline
                'fprintf(''       TRANSIENT STABILITY RESULTS\n'')' newline
                'fprintf(''========================================\n'')' newline
                'fprintf(''Initial rotor angle: %.2f degrees\n'', delta0*180/pi);' newline
                'fprintf(''Maximum rotor angle: %.2f degrees\n'', max(delta)*180/pi);' newline
                'fprintf(''Final rotor angle: %.2f degrees\n'', delta(end)*180/pi);' newline
                'fprintf(''Maximum speed deviation: %.4f pu\n'', max(abs(omega)));' newline
                newline
                '% Stability assessment' newline
                'if max(delta) < pi/2' newline
                '    fprintf(''\nSystem is STABLE\n'');' newline
                'else' newline
                '    fprintf(''\nSystem is UNSTABLE - Loss of synchronism\n'');' newline
                'end' newline
                'fprintf(''========================================\n'')' newline
                newline
                '%% Rotor Angle Swing Curve' newline
                'figure;' newline
                'plot(t, delta*180/pi, ''b-'', ''LineWidth'', 2);' newline
                'hold on;' newline
                'xline(tFault, ''r--'', ''Fault Applied'');' newline
                'xline(tClear, ''g--'', ''Fault Cleared'');' newline
                'yline(90, ''k:'', ''Critical Angle'');' newline
                'xlabel(''Time (seconds)'');' newline
                'ylabel(''Rotor Angle (degrees)'');' newline
                'title(''Rotor Angle Swing Curve'');' newline
                'grid on;' newline
                newline
                '%% Speed Deviation' newline
                'figure;' newline
                'plot(t, omega, ''b-'', ''LineWidth'', 2);' newline
                'xlabel(''Time (seconds)'');' newline
                'ylabel(''Speed Deviation (pu)'');' newline
                'title(''Generator Speed Deviation'');' newline
                'grid on;' newline
                newline
                '%% Phase Plane Trajectory' newline
                'figure;' newline
                'plot(delta*180/pi, omega, ''b-'', ''LineWidth'', 1.5);' newline
                'hold on;' newline
                'plot(delta(1)*180/pi, omega(1), ''go'', ''MarkerFaceColor'', ''g'');' newline
                'plot(delta(end)*180/pi, omega(end), ''ro'', ''MarkerFaceColor'', ''r'');' newline
                'xlabel(''Rotor Angle (degrees)'');' newline
                'ylabel(''Speed Deviation (pu)'');' newline
                'title(''Phase Plane Trajectory'');' newline
                'grid on;' newline
                newline
            ];

            % Add critical clearing time calculation
            script = [script
                '%% Critical Clearing Time Study' newline
                'fprintf(''\n--- Critical Clearing Time Analysis ---\n'');' newline
                'clearingTimes = 0.05:0.01:0.3;' newline
                'isStable = zeros(size(clearingTimes));' newline
                newline
                'for k = 1:length(clearingTimes)' newline
                '    tClear_test = clearingTimes(k);' newline
                '    delta_test = zeros(1, nSteps);' newline
                '    omega_test = zeros(1, nSteps);' newline
                '    delta_test(1) = delta0;' newline
                '    omega_test(1) = omega0;' newline
                '    ' newline
                '    for i = 2:nSteps' newline
                '        if t(i) < tFault' newline
                '            Pmax_current = Pmax_prefault;' newline
                '        elseif t(i) < tClear_test' newline
                '            Pmax_current = Pmax_during_fault;' newline
                '        else' newline
                '            Pmax_current = Pmax_post_fault;' newline
                '        end' newline
                '        ' newline
                '        Pe = Pmax_current * sin(delta_test(i-1));' newline
                '        domega_dt = (Pm - Pe - D * omega_test(i-1)) / M;' newline
                '        ddelta_dt = omega_test(i-1);' newline
                '        ' newline
                '        omega_test(i) = omega_test(i-1) + domega_dt * dt;' newline
                '        delta_test(i) = delta_test(i-1) + ddelta_dt * dt;' newline
                '    end' newline
                '    ' newline
                '    isStable(k) = max(delta_test) < pi/2;' newline
                'end' newline
                newline
                'criticalIdx = find(isStable == 0, 1);' newline
                'if ~isempty(criticalIdx)' newline
                '    CCT = clearingTimes(criticalIdx - 1);' newline
                '    fprintf(''Critical Clearing Time: %.3f seconds\n'', CCT);' newline
                'else' newline
                '    fprintf(''System stable for all tested clearing times\n'');' newline
                '    CCT = max(clearingTimes);' newline
                'end' newline
                newline
                'figure;' newline
                'plot(clearingTimes, isStable, ''bo-'', ''LineWidth'', 2);' newline
                'xlabel(''Clearing Time (seconds)'');' newline
                'ylabel(''Stable (1) / Unstable (0)'');' newline
                'title(''Critical Clearing Time Analysis'');' newline
                'grid on;' newline
            ];
        end

        function script = generateSmallSignalCode(obj)
            % Generate small signal stability analysis code

            script = [
                '%% Small Signal Stability Analysis' newline
                '% Auto-generated MATLAB code for eigenvalue analysis' newline
                newline
                'clear; clc; close all;' newline
                newline
                '%% System Parameters' newline
                sprintf('f = 60;  %% System frequency (Hz)' newline)
                sprintf('omega_s = 2*pi*f;  %% Synchronous speed (rad/s)' newline)
                newline
                '%% Generator Parameters' newline
            ];

            if isempty(obj.generatorData)
                script = [script
                    'H = 4.0;  %% Inertia constant (s)' newline
                    'D = 0.01;  %% Damping coefficient (pu)' newline
                    'K1 = 0.8;  %% Synchronizing torque coefficient' newline
                    'K2 = 1.2;  %% Damping torque coefficient' newline
                ];
            else
                gen = obj.generatorData(1);
                script = [script
                    sprintf('H = %g;  %% Inertia constant (s)\n', gen.H)
                    sprintf('D = %g;  %% Damping coefficient (pu)\n', gen.D)
                    'K1 = 0.8;  %% Synchronizing torque coefficient' newline
                    'K2 = 1.2;  %% Damping torque coefficient' newline
                ];
            end

            script = [script
                newline
                '%% State Matrix Construction' newline
                '% Linearized swing equation:' newline
                '% d(delta)/dt = omega' newline
                '% d(omega)/dt = (Pm - Pe - D*omega)/M' newline
                '% Linearized: d(Delta_omega)/dt = -K1/M * Delta_delta - D/M * Delta_omega' newline
                newline
                'M = 2 * H / omega_s;  %% Inertia coefficient' newline
                newline
                '% State matrix A' newline
                'A = [0, 1;' newline
                '     -K1/M, -D/M];' newline
                newline
                '%% Eigenvalue Analysis' newline
                '[V, D_eig] = eig(A);' newline
                'eigenvalues = diag(D_eig);' newline
                newline
                'fprintf(''\n========================================\n'')' newline
                'fprintf(''    SMALL SIGNAL STABILITY RESULTS\n'')' newline
                'fprintf(''========================================\n'')' newline
                'fprintf(''\nEigenvalues:\n'');' newline
                'for i = 1:length(eigenvalues)' newline
                '    fprintf(''  Lambda_%d = %.4f + j%.4f\n'', i, real(eigenvalues(i)), imag(eigenvalues(i)));' newline
                'end' newline
                newline
                '%% Calculate Modal Parameters' newline
                'for i = 1:length(eigenvalues)' newline
                '    lambda = eigenvalues(i);' newline
                '    if imag(lambda) ~= 0' newline
                '        omega_n(i) = abs(lambda);  %% Natural frequency (rad/s)' newline
                '        f_n(i) = omega_n(i) / (2*pi);  %% Frequency (Hz)' newline
                '        zeta(i) = -real(lambda) / omega_n(i);  %% Damping ratio' newline
                '        sigma(i) = real(lambda);  %% Decay rate' newline
                '    else' newline
                '        omega_n(i) = abs(lambda);' newline
                '        f_n(i) = 0;' newline
                '        zeta(i) = 1;  %% Overdamped' newline
                '        sigma(i) = real(lambda);' newline
                '    end' newline
                'end' newline
                newline
                'fprintf(''\n--- Modal Analysis ---\n'');' newline
                'for i = 1:length(eigenvalues)' newline
                '    if imag(eigenvalues(i)) ~= 0' newline
                '        fprintf(''Mode %d:\n'', i);' newline
                '        fprintf(''  Natural frequency: %.2f Hz\n'', f_n(i));' newline
                '        fprintf(''  Damping ratio: %.2f %%\n'', zeta(i)*100);' newline
                '        if zeta(i) < 0.05' newline
                '            fprintf(''  WARNING: Poorly damped mode!\n'');' newline
                '        end' newline
                '    end' newline
                'end' newline
                newline
                '%% Stability Assessment' newline
                'if all(real(eigenvalues) < 0)' newline
                '    fprintf(''\nSystem is SMALL SIGNAL STABLE\n'');' newline
                'else' newline
                '    fprintf(''\nSystem is SMALL SIGNAL UNSTABLE\n'');' newline
                'end' newline
                'fprintf(''========================================\n'')' newline
                newline
                '%% Eigenvalue Plot' newline
                'figure;' newline
                'plot(real(eigenvalues), imag(eigenvalues), ''bx'', ''MarkerSize'', 10, ''LineWidth'', 2);' newline
                'hold on;' newline
                'plot(real(eigenvalues), -imag(eigenvalues), ''bx'', ''MarkerSize'', 10, ''LineWidth'', 2);' newline
                'xline(0, ''r--'', ''Stability Boundary'');' newline
                'xlabel(''Real (1/s)'');' newline
                'ylabel(''Imaginary (rad/s)'');' newline
                'title(''Eigenvalue Plot'');' newline
                'grid on;' newline
                'axis equal;' newline
                newline
                '%% Damping Ratio vs Frequency' newline
                'figure;' newline
                'plot(f_n, zeta*100, ''bo'', ''MarkerFaceColor'', ''b'', ''MarkerSize'', 8);' newline
                'yline(5, ''r--'', ''5%% Damping'');' newline
                'xlabel(''Frequency (Hz)'');' newline
                'ylabel(''Damping Ratio (%%)'');' newline
                'title(''Damping Ratio vs Oscillation Frequency'');' newline
                'grid on;' newline
            ];
        end

        function script = generateVoltageStabilityCode(obj)
            % Generate voltage stability analysis code

            script = [
                '%% Voltage Stability Analysis' newline
                '% Auto-generated MATLAB code for P-V and Q-V curve analysis' newline
                newline
                'clear; clc; close all;' newline
                newline
                '%% System Parameters' newline
                'V_source = 1.0;  %% Source voltage (pu)' newline
                'X = 0.5;  %% System reactance (pu)' newline
                'R = 0.05;  %% System resistance (pu)' newline
                'Z = R + 1i*X;  %% System impedance' newline
                newline
                '%% P-V Curve Analysis' newline
                '% For a simple two-bus system:' newline
                '% P = V^2/R * (1 - V/V_source) for constant power factor' newline
                newline
                'V_range = linspace(0, 1.2, 500);  %% Voltage range' newline
                'pf = 0.95;  %% Power factor' newline
                'theta = acos(pf);' newline
                newline
                '% Maximum power transfer' newline
                'Pmax = V_source^2 / (4 * X);  %% Simplified for lossless case' newline
                newline
                '% Calculate P-V relationship' newline
                'P_PV = zeros(size(V_range));' newline
                'for i = 1:length(V_range)' newline
                '    V = V_range(i);' newline
                '    % Power flow equation for two-bus system' newline
                '    delta = acos((V.^2 + V_source^2 - 2*X*P_PV(i))./(2*V*V_source));' newline
                'end' newline
                newline
                '% Simplified P-V curve (nose curve)' newline
                'P_PV = (V_source * V_range ./ X) .* sin(acos(V_range / V_source));' newline
                'P_PV(imag(P_PV) ~= 0) = NaN;  %% Remove invalid points' newline
                newline
                '%% Q-V Curve Analysis' newline
                'Q_QV = (V_source * V_range ./ X) .* cos(acos(V_range / V_source)) - V_range.^2 / X;' newline
                'Q_QV(imag(Q_QV) ~= 0) = NaN;' newline
                newline
                'fprintf(''\n========================================\n'')' newline
                'fprintf(''      VOLTAGE STABILITY RESULTS\n'')' newline
                'fprintf(''========================================\n'')' newline
                'fprintf(''Maximum Power Transfer: %.3f pu\n'', Pmax);' newline
                'fprintf(''Critical Voltage: %.3f pu\n'', 0.707);' newline
                'fprintf(''Voltage Stability Margin: %.1f %%\n'', (1 - 0.707) * 100);' newline
                'fprintf(''========================================\n'')' newline
                newline
                '%% P-V Curve Plot' newline
                'figure;' newline
                'subplot(2,1,1);' newline
                'plot(real(P_PV), V_range, ''b-'', ''LineWidth'', 2);' newline
                'hold on;' newline
                'plot(Pmax, 0.707, ''ro'', ''MarkerFaceColor'', ''r'', ''MarkerSize'', 8);' newline
                'xlabel(''Active Power P (pu)'');' newline
                'ylabel(''Voltage Magnitude (pu)'');' newline
                'title(''P-V Curve (Nose Curve)'');' newline
                'grid on;' newline
                'text(Pmax, 0.75, ''Critical Point'', ''HorizontalAlignment'', ''center'');' newline
                newline
                '%% Q-V Curve Plot' newline
                'subplot(2,1,2);' newline
                'plot(real(Q_QV), V_range, ''r-'', ''LineWidth'', 2);' newline
                'xlabel(''Reactive Power Q (pu)'');' newline
                'ylabel(''Voltage Magnitude (pu)'');' newline
                'title(''Q-V Curve'');' newline
                'grid on;' newline
                newline
                '%% Voltage Stability Indices' newline
                '% L-index calculation' newline
                'V_operating = 1.0;  %% Normal operating voltage' newline
                'V_critical = 0.707;  %% Critical voltage' newline
                'L_index = 1 - (V_operating / V_critical)^2;' newline
                newline
                'fprintf(''\n--- Voltage Stability Indices ---\n'');' newline
                'fprintf(''L-Index: %.3f\n'', L_index);' newline
                'if L_index < 0.2' newline
                '    fprintf(''Status: Good voltage stability margin\n'');' newline
                'elseif L_index < 0.4' newline
                '    fprintf(''Status: Adequate margin\n'');' newline
                'else' newline
                '    fprintf(''Status: Warning - Low stability margin\n'');' newline
                'end' newline
            ];
        end

        function script = generateFrequencyStabilityCode(obj)
            % Generate frequency stability analysis code

            script = [
                '%% Frequency Stability Analysis' newline
                '% Auto-generated MATLAB code for load-frequency dynamics' newline
                newline
                'clear; clc; close all;' newline
                newline
                '%% System Parameters' newline
                sprintf('f0 = 60;  %% Nominal frequency (Hz)' newline)
                sprintf('simTime = %g;  %% Simulation time (s)' newline, obj.simulationTime)
                newline
                '%% Generator Parameters' newline
                'H = 4.0;  %% Inertia constant (s)' newline
                'R = 0.05;  %% Speed regulation (pu)' newline
                'D = 0.01;  %% Load damping (pu)' newline
                'Tg = 0.2;  %% Governor time constant (s)' newline
                'Tt = 0.5;  %% Turbine time constant (s)' newline
                newline
                '%% Load Disturbance' newline
                'dPL = 0.1;  %% Load step increase (pu)' newline
                newline
                '%% Simulation' newline
                'dt = 0.01;  %% Time step' newline
                't = 0:dt:simTime;' newline
                'nSteps = length(t);' newline
                newline
                '% State variables' newline
                'df = zeros(1, nSteps);  %% Frequency deviation' newline
                'dPg = zeros(1, nSteps);  %% Generator power change' newline
                'dx = zeros(1, nSteps);  %% Governor valve position' newline
                newline
                '%% Simulation Loop' newline
                'for i = 2:nSteps' newline
                '    % Load disturbance (step at t=1s)' newline
                '    if t(i) >= 1' newline
                '        dPL_current = dPL;' newline
                '    else' newline
                '        dPL_current = 0;' newline
                '    end' newline
                '    ' newline
                '    % Governor dynamics' newline
                '    dx_dt = (-df(i-1)/R - dx(i-1)) / Tg;' newline
                '    ' newline
                '    % Turbine dynamics' newline
                '    dPg_dt = (dx(i-1) - dPg(i-1)) / Tt;' newline
                '    ' newline
                '    % Swing equation for frequency' newline
                '    df_dt = (dPg(i-1) - dPL_current - D * df(i-1)) / (2*H);' newline
                '    ' newline
                '    % Euler integration' newline
                '    dx(i) = dx(i-1) + dx_dt * dt;' newline
                '    dPg(i) = dPg(i-1) + dPg_dt * dt;' newline
                '    df(i) = df(i-1) + df_dt * dt;' newline
                'end' newline
                newline
                '%% Display Results' newline
                'fprintf(''\n========================================\n'')' newline
                'fprintf(''    FREQUENCY STABILITY RESULTS\n'')' newline
                'fprintf(''========================================\n'')' newline
                'fprintf(''Load disturbance: %.2f pu\n'', dPL);' newline
                'fprintf(''Maximum frequency deviation: %.3f Hz\n'', min(df)*f0);' newline
                'fprintf(''Steady-state frequency deviation: %.3f Hz\n'', df(end)*f0);' newline
                'fprintf(''Minimum frequency: %.2f Hz\n'', f0 + min(df)*f0);' newline
                'fprintf(''========================================\n'')' newline
                newline
                '%% Frequency Response Plot' newline
                'figure;' newline
                'plot(t, f0 + df*f0, ''b-'', ''LineWidth'', 2);' newline
                'xlabel(''Time (seconds)'');' newline
                'ylabel(''Frequency (Hz)'');' newline
                'title(''Frequency Response to Load Change'');' newline
                'grid on;' newline
                'yline(f0, ''k--'', ''Nominal'');' newline
                'yline(59.5, ''r--'', ''Lower Limit'');' newline
                newline
                '%% Generator Power Response' newline
                'figure;' newline
                'plot(t, dPg, ''g-'', ''LineWidth'', 2);' newline
                'xlabel(''Time (seconds)'');' newline
                'ylabel(''Generator Power Change (pu)'');' newline
                'title(''Governor-Turbine Response'');' newline
                'grid on;' newline
            ];
        end

        function generateSimulinkModel(obj, params)
            % Generate Simulink model for stability simulation
            fprintf('\nGenerating Simulink stability model...\n\n');

            script = [
                '% Transient Stability Simulation Model' newline
                'modelName = ''transient_stability'';' newline
                'if system_exists(modelName)' newline
                '    close_system(modelName, 0);' newline
                '    delete_system(modelName);' newline
                'end' newline
                'new_system(modelName);' newline
                newline
                '% Add synchronous machine block' newline
                'add_block(''powerlib/Machines/Simplified Synchronous Machine'', ...' newline
                '    [modelName ''/Synchronous Machine'']);' newline
                newline
                '% Add governor and exciter' newline
                'add_block(''powerlib/Control Blocks/PI Controller'', ...' newline
                '    [modelName ''/Governor'']);' newline
                newline
                '% Add fault block' newline
                'add_block(''powerlib/Faults/Three-Phase Fault'', ...' newline
                '    [modelName ''/Fault'']);' newline
                newline
                '% Add scopes' newline
                'add_block(''simulink/Sinks/Scope'', [modelName ''/Rotor Angle'']);' newline
                'add_block(''simulink/Sinks/Scope'', [modelName ''/Speed'']);' newline
                newline
                'open_system(modelName);' newline
                'fprintf(''Stability simulation model created: %s.slx\n'', modelName);' newline
            ];

            fprintf('%s\n', script);
        end

        function handleInput(obj, userInput)
            % Handle contextual user input
            fprintf('Stability module received: %s\n', userInput);
        end
    end
end
