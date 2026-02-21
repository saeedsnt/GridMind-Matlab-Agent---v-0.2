classdef FaultAnalysisModule < handle
    % FaultAnalysisModule - Fault analysis domain module
    %
    % Provides interactive parameter collection and code generation
    % for short circuit studies and fault calculations.

    properties
        faultType             % Type of fault
        faultLocation         % Bus number for fault
        faultImpedance        % Fault impedance (pu)
        XR_ratio              % X/R ratio of system
        voltageLevel          % System voltage level (kV)
        baseMVA               % Base MVA
        prefaultVoltage       % Pre-fault voltage (pu)
        includeDCoffset       % Include DC offset in calculations
        busData               % Bus data
        lineData              % Line data
        generatorData         % Generator data
    end

    methods
        function obj = FaultAnalysisModule()
            % Constructor - set defaults
            obj.faultType = 'three_phase';
            obj.faultLocation = 1;
            obj.faultImpedance = 0.01;
            obj.XR_ratio = 10;
            obj.voltageLevel = 138;
            obj.baseMVA = 100;
            obj.prefaultVoltage = 1.0;
            obj.includeDCoffset = true;
            obj.busData = [];
            obj.lineData = [];
            obj.generatorData = [];
        end

        function params = collectParameters(obj, prompt)
            % Collect fault analysis parameters from user

            fprintf('\n=== Fault Analysis Configuration ===\n\n');

            % Select fault type
            faultTypes = {
                'three_phase', 'single_line_ground', ...
                'line_to_line', 'double_line_ground'
            };
            faultDescriptions = {
                'Three-phase fault (symmetrical)', ...
                'Single line-to-ground fault', ...
                'Line-to-line fault', ...
                'Double line-to-ground fault'
            };

            fprintf('Available fault types:\n');
            faultIdx = prompt.getSelectOne(faultTypes, faultDescriptions, 'Select fault type:');
            obj.faultType = faultIdx;
            fprintf('Selected fault type: %s\n\n', obj.faultType);

            % System parameters
            obj.voltageLevel = prompt.getNumericInputWithDefault( ...
                'Enter system voltage level (kV): ', 138, [0.12, 1000]);

            obj.baseMVA = prompt.getNumericInputWithDefault( ...
                'Enter base MVA: ', 100, [1, 10000]);

            % Fault location
            obj.numBuses = prompt.getNumericInput('Enter number of buses in system: ', [2, 500]);
            obj.faultLocation = prompt.getNumericInput( ...
                'Enter fault location (bus number): ', [1, obj.numBuses]);

            % Fault parameters
            obj.faultImpedance = prompt.getNumericInputWithDefault( ...
                'Enter fault impedance (pu): ', 0.01, [0, 1]);

            % Determine typical X/R ratio based on voltage
            if obj.voltageLevel >= 230
                defaultXR = 15;
            elseif obj.voltageLevel >= 69
                defaultXR = 12;
            elseif obj.voltageLevel >= 34.5
                defaultXR = 8;
            else
                defaultXR = 3;
            end
            obj.XR_ratio = prompt.getNumericInputWithDefault( ...
                'Enter X/R ratio: ', defaultXR, [0.1, 50]);

            % Pre-fault conditions
            obj.prefaultVoltage = prompt.getNumericInputWithDefault( ...
                'Pre-fault voltage at fault location (pu): ', 1.0, [0.9, 1.1]);

            % DC offset
            obj.includeDCoffset = prompt.getYesNo( ...
                'Include DC offset in fault current calculation? (y/n): ');

            % Collect system data
            obj.busData = obj.collectBusData(prompt, obj.numBuses);
            obj.lineData = obj.collectLineData(prompt, obj.numBuses);
            obj.generatorData = obj.collectGeneratorData(prompt);

            % Return parameters struct
            params.faultType = obj.faultType;
            params.faultLocation = obj.faultLocation;
            params.faultImpedance = obj.faultImpedance;
            params.XR_ratio = obj.XR_ratio;
            params.voltageLevel = obj.voltageLevel;
            params.baseMVA = obj.baseMVA;
            params.prefaultVoltage = obj.prefaultVoltage;
            params.includeDCoffset = obj.includeDCoffset;
            params.numBuses = obj.numBuses;
            params.busData = obj.busData;
            params.lineData = obj.lineData;
            params.generatorData = obj.generatorData;
        end

        function busData = collectBusData(obj, prompt, numBuses)
            % Collect simplified bus data for fault analysis
            busData = struct('Z1', {}, 'Z0', {}, 'type', {});

            fprintf('\n--- Bus Impedance Data ---\n');
            fprintf('Enter positive and zero sequence impedances\n\n');

            for i = 1:numBuses
                fprintf('Bus %d:\n', i);
                busData(i).Z1 = prompt.getNumericInputWithDefault( ...
                    '  Positive sequence impedance Z1 (pu): ', 0.01 + 0.05i, [0, 10]);
                busData(i).Z0 = prompt.getNumericInputWithDefault( ...
                    '  Zero sequence impedance Z0 (pu): ', 0.03 + 0.15i, [0, 10]);
                busData(i).type = 'bus';
                fprintf('\n');
            end
        end

        function lineData = collectLineData(obj, prompt, numBuses)
            % Collect line data for fault analysis
            lineData = struct('fromBus', {}, 'toBus', {}, 'Z1', {}, 'Z0', {});

            fprintf('\n--- Line Sequence Impedance Data ---\n');
            numLines = prompt.getNumericInput('Enter number of lines: ', [1, numBuses*(numBuses-1)]);

            for i = 1:numLines
                fprintf('Line %d:\n', i);
                lineData(i).fromBus = prompt.getNumericInput('  From bus: ', [1, numBuses]);
                lineData(i).toBus = prompt.getNumericInput('  To bus: ', [1, numBuses]);
                lineData(i).Z1 = prompt.getNumericInputWithDefault( ...
                    '  Positive sequence impedance Z1 (pu): ', 0.02 + 0.08i, [0, 10]);
                lineData(i).Z0 = prompt.getNumericInputWithDefault( ...
                    '  Zero sequence impedance Z0 (pu): ', 0.06 + 0.24i, [0, 10]);
                fprintf('\n');
            end
        end

        function genData = collectGeneratorData(obj, prompt)
            % Collect generator data for fault analysis
            genData = struct('bus', {}, 'Xd', {}, 'Xdp', {}, 'Xdp2', {}, 'Ra', {});

            fprintf('\n--- Generator Data ---\n');
            numGens = prompt.getNumericInput('Enter number of generators: ', [1, 50]);

            for i = 1:numGens
                fprintf('Generator %d:\n', i);
                genData(i).bus = prompt.getNumericInput('  Connected bus: ', [1, obj.numBuses]);
                genData(i).Xd = prompt.getNumericInputWithDefault( ...
                    '  Synchronous reactance Xd (pu): ', 1.8, [0.1, 10]);
                genData(i).Xdp = prompt.getNumericInputWithDefault( ...
                    '  Transient reactance X''d (pu): ', 0.3, [0.05, 2]);
                genData(i).Xdp2 = prompt.getNumericInputWithDefault( ...
                    '  Subtransient reactance X''''d (pu): ', 0.2, [0.05, 1]);
                genData(i).Ra = prompt.getNumericInputWithDefault( ...
                    '  Armature resistance Ra (pu): ', 0.003, [0, 0.1]);
                fprintf('\n');
            end
        end

        function code = generateCode(obj, params)
            % Generate MATLAB code for fault analysis
            code = obj.generateFaultScript();
        end

        function script = generateFaultScript(obj)
            % Generate complete fault analysis script

            script = [
                '%% Fault Analysis - Short Circuit Study' newline
                '% Auto-generated MATLAB code for fault calculations' newline
                newline
                'clear; clc; close all;' newline
                newline
                '%% System Parameters' newline
                sprintf('baseMVA = %g;' newline, obj.baseMVA)
                sprintf('voltageLevel = %g;  %% kV' newline, obj.voltageLevel)
                sprintf('faultType = '''', obj.faultType)
                ''';' newline
                sprintf('faultLocation = %d;' newline, obj.faultLocation)
                sprintf('Zf = %g;  %% Fault impedance (pu)' newline, obj.faultImpedance)
                sprintf('XR_ratio = %g;' newline, obj.XR_ratio)
                sprintf('Vprefault = %g;  %% Pre-fault voltage (pu)' newline, obj.prefaultVoltage)
                newline
                '%% Base Values' newline
                'baseKV = voltageLevel;' newline
                'baseZ = baseKV^2 / baseMVA;' newline
                'baseI = baseMVA * 1000 / (sqrt(3) * baseKV);  %% A' newline
                newline
            ];

            % Add generator data
            script = [script
                '%% Generator Data (per-unit on system base)' newline
                '% Format: [Bus, Xd, X''d, X''''d, Ra]' newline
                'genData = [' newline
            ];

            if isempty(obj.generatorData)
                % Default generator at bus 1
                script = [script, '    1, 1.8, 0.3, 0.2, 0.003;' newline];
            else
                for i = 1:length(obj.generatorData)
                    gen = obj.generatorData(i);
                    script = [script, sprintf('    %d, %.3f, %.3f, %.3f, %.4f', ...
                        gen.bus, gen.Xd, gen.Xdp, gen.Xdp2, gen.Ra)];
                    if i < length(obj.generatorData)
                        script = [script, ';' newline];
                    else
                        script = [script, ';' newline];
                    end
                end
            end
            script = [script, '];' newline newline];

            % Add network impedance data
            script = [script
                '%% Network Impedance Data' newline
                '% Positive sequence impedance matrix (simplified)' newline
                'Z1 = [' newline
            ];

            % Build impedance matrix representation
            if isempty(obj.busData)
                script = [script, '    0.01 + 0.05i, 0.005 + 0.02i, 0.005 + 0.02i;' newline];
                script = [script, '    0.005 + 0.02i, 0.01 + 0.05i, 0.005 + 0.02i;' newline];
                script = [script, '    0.005 + 0.02i, 0.005 + 0.02i, 0.01 + 0.05i' newline];
            else
                for i = 1:min(length(obj.busData), 5)
                    bus = obj.busData(i);
                    script = [script, sprintf('    %s', mat2str(bus.Z1))];
                    if i < min(length(obj.busData), 5)
                        script = [script, ';' newline];
                    end
                end
            end
            script = [script, '];' newline newline];

            script = [script
                '% Zero sequence impedance (typically 2-3 times positive sequence)' newline
                'Z0 = 3 * Z1;' newline
                newline
            ];

            % Add fault calculation code based on fault type
            script = [script, obj.getFaultCalculationsCode()];

            % Add results display
            script = [script
                newline
                '%% Display Results' newline
                'fprintf(''\n========================================\n'')' newline
                'fprintf(''       FAULT ANALYSIS RESULTS\n'')' newline
                'fprintf(''========================================\n'')' newline
                'fprintf(''Fault Type: %s\n'', faultType)' newline
                'fprintf(''Fault Location: Bus %d\n'', faultLocation)' newline
                'fprintf(''Fault Impedance: %.3f pu\n'', Zf)' newline
                'fprintf(''\n--- Fault Currents ---\n'')' newline
                sprintf('fprintf(''Fault Current: %.3f pu (%.2f kA)\n'', Ifault_pu, Ifault_kA);' newline)
                'fprintf(''DC Component (initial): %.2f kA\n'', Idc_initial);' newline
                'fprintf(''Total Asymmetrical: %.2f kA\n'', Ifault_asym);' newline
                newline
                'fprintf(''\n--- Breaker Rating Requirements ---\n'')' newline
                'fprintf(''Interrupting Capacity: %.2f kA\n'', Iinterrupting);' newline
                'fprintf(''Momentary (Close & Latch): %.2f kA\n'', Imomentary);' newline
                'fprintf(''\n========================================\n'')' newline
            ];

            % Add plots
            script = [script
                newline
                '%% Fault Current Waveform' newline
                'figure;' newline
                'plot(t*1000, ifault_ac + ifault_dc, ''b'', ''LineWidth'', 1.5);' newline
                'hold on;' newline
                'plot(t*1000, ifault_ac, ''r--'', ''LineWidth'', 1);' newline
                'plot(t*1000, -ifault_ac, ''r--'', ''LineWidth'', 1);' newline
                'xlabel(''Time (ms)'');' newline
                'ylabel(''Current (kA)'');' newline
                'title([''Fault Current Waveform - '' faultType]);' newline
                'legend(''Total Current'', ''AC Component'');' newline
                'grid on;' newline
                newline
                '%% DC Offset Decay' newline
                'figure;' newline
                'plot(t*1000, ifault_dc, ''b'', ''LineWidth'', 1.5);' newline
                'xlabel(''Time (ms)'');' newline
                'ylabel(''DC Component (kA)'');' newline
                'title(''DC Offset Decay'');' newline
                'grid on;' newline
            ];
        end

        function code = getFaultCalculationsCode(obj)
            % Generate fault calculation code based on fault type

            code = [
                '%% Fault Current Calculations' newline
                newline
                '% System parameters' newline
                'RoverX = 1 / XR_ratio;' newline
                'X1 = imag(Z1(faultLocation, faultLocation));' newline
                'R1 = R1 = X1 * RoverX;' newline
                'Z1_eq = R1 + 1i*X1;' newline
                newline
                '% Time constants' newline
                'f = 60;  %% Hz' newline
                'omega = 2*pi*f;' newline
                'tau_dc = XR_ratio / omega;  %% DC time constant' newline
                newline
                '% Symmetrical fault current calculation based on fault type' newline
                'switch faultType' newline
                '    case ''three_phase''' newline
                '        %% Three-phase fault (balanced)' newline
                '        If_sym = Vprefault / (Z1_eq + Zf);' newline
                '        Ifault_pu = abs(If_sym);' newline
                '        fprintf(''\nThree-phase fault - maximum fault current\n'');' newline
                newline
                '    case ''single_line_ground''' newline
                '        %% Single line-to-ground fault' newline
                '        %% I_f = 3*V / (Z1 + Z2 + Z0 + 3*Zf)' newline
                '        Z2 = Z1_eq;  %% Assume Z2 = Z1' newline
                '        If_sym = 3 * Vprefault / (Z1_eq + Z2 + Z0(faultLocation,faultLocation) + 3*Zf);' newline
                '        Ifault_pu = abs(If_sym);' newline
                '        fprintf(''\nSLG fault - most common fault type\n'');' newline
                newline
                '    case ''line_to_line''' newline
                '        %% Line-to-line fault' newline
                '        %% I_f = sqrt(3)*V / (Z1 + Z2 + Zf)' newline
                '        Z2 = Z1_eq;' newline
                '        If_sym = sqrt(3) * Vprefault / (Z1_eq + Z2 + Zf);' newline
                '        Ifault_pu = abs(If_sym);' newline
                '        fprintf(''\nLine-to-line fault\n'');' newline
                newline
                '    case ''double_line_ground''' newline
                '        %% Double line-to-ground fault' newline
                '        Z2 = Z1_eq;' newline
                '        Z0_net = Z0(faultLocation,faultLocation);' newline
                '        %% I_f = sqrt(3)*V * sqrt(1 - (Z2*Z0)/((Z2+Z0)^2)) / (Z1 + Z2*Z0/(Z2+Z0) + Zf)' newline
                '        Z_parallel = (Z2 * Z0_net) / (Z2 + Z0_net);' newline
                '        If_sym = sqrt(3) * Vprefault / (Z1_eq + Z_parallel + Zf);' newline
                '        Ifault_pu = abs(If_sym);' newline
                '        fprintf(''\nDouble line-to-ground fault\n'');' newline
                'end' newline
                newline
                '% Convert to actual values' newline
                'Ifault_kA = Ifault_pu * baseI / 1000;  %% kA' newline
                newline
                '%% DC Offset Calculation' newline
                '% Maximum DC offset occurs when fault at voltage zero crossing' newline
                'Idc_initial = sqrt(2) * Ifault_kA;  %% Initial DC component' newline
                't = linspace(0, 0.5, 1000);  %% 500 ms' newline
                'ifault_dc = Idc_initial * exp(-t / tau_dc);  %% kA' newline
                newline
                '%% AC Component with Decay' newline
                '% Simplified decay model' newline
                'Ifault_asym = sqrt(Ifault_kA^2 + Idc_initial^2);  %% Initial asymmetrical' newline
                'ifault_ac = sqrt(2) * Ifault_kA * sin(omega*t);  %% kA' newline
                newline
                '%% Total Fault Current' newline
                'ifault_total = ifault_ac + ifault_dc;  %% kA' newline
                newline
                '%% Breaker Rating Calculations' newline
                '% Interrupting capacity (symmetrical)' newline
                'Iinterrupting = Ifault_kA;  %% kA' newline
                newline
                '% Momentary rating (close & latch) - includes DC offset' newline
                '% Typically 1.6 to 2.7 times symmetrical depending on X/R' newline
                'asymmetryFactor = 1 + exp(-0.5/(3*tau_dc));  %% At 0.5 cycles' newline
                'Imomentary = Ifault_kA * asymmetryFactor * sqrt(2);  %% Peak close & latch' newline
                newline
                '%% Contact Parting Current (for breaker selection)' newline
                '% Assuming 3-cycle (50ms) breaker' newline
                't_contact = 0.05;  %% 3 cycles at 60 Hz' newline
                'Idc_at_contact = Idc_initial * exp(-t_contact / tau_dc);' newline
                'Icontact = sqrt(Ifault_kA^2 + (Idc_at_contact/sqrt(2))^2);' newline
                'fprintf(''Contact Parting Current: %.2f kA\n'', Icontact);' newline
            ];
        end

        function generateSimulinkModel(obj, params)
            % Generate Simulink model for fault simulation
            fprintf('\nGenerating Simulink fault simulation model...\n\n');

            script = [
                '% Fault Simulation Model' newline
                'modelName = ''fault_simulation'';' newline
                'if system_exists(modelName)' newline
                '    close_system(modelName, 0);' newline
                '    delete_system(modelName);' newline
                'end' newline
                'new_system(modelName);' newline
                newline
                '% Add three-phase source' newline
                'add_block(''powerlib/Electrical Sources/Three-Phase Source'', ...' newline
                '    [modelName ''/Three-Phase Source'']);' newline
                newline
                '% Add three-phase fault block' newline
                'add_block(''powerlib/Faults/Three-Phase Fault'', ...' newline
                '    [modelName ''/Three-Phase Fault'']);' newline
                newline
                '% Add measurement blocks' newline
                'add_block(''powerlib/Measurements/V-I Measurement'', ...' newline
                '    [modelName ''/VI Measurement'']);' newline
                newline
                '% Connect blocks and set parameters' newline
                '% (Full implementation requires Simscape Electrical)' newline
                newline
                'open_system(modelName);' newline
                'fprintf(''Fault simulation model created: %s.slx\n'', modelName);' newline
            ];

            fprintf('%s\n', script);
        end

        function handleInput(obj, userInput)
            % Handle contextual user input
            fprintf('Fault analysis module received: %s\n', userInput);
        end
    end
end
