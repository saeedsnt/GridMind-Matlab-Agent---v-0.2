classdef ProtectionModule < handle
    % ProtectionModule - Protection systems domain module
    %
    % Provides interactive parameter collection and code generation
    % for relay coordination and protection system design.

    properties
        relayType             % Type of relay
        CT_ratio              % Current transformer ratio
        VT_ratio              % Voltage transformer ratio
        pickupCurrent         % Relay pickup current
        TDS                   % Time dial setting
        CTI                   % Coordination time interval
        curveType             % Relay curve type
        systemType            % Distribution or transmission
        zones                 % Protection zones for distance relays
    end

    methods
        function obj = ProtectionModule()
            % Constructor - set defaults
            obj.relayType = 'overcurrent';
            obj.CT_ratio = 400;
            obj.VT_ratio = 120;
            obj.pickupCurrent = 100;
            obj.TDS = 1.0;
            obj.CTI = 0.3;
            obj.curveType = 'IEEE moderately inverse';
            obj.systemType = 'distribution';
            obj.zones = struct('reach', {}, 'time', {});
        end

        function params = collectParameters(obj, prompt)
            % Collect protection system parameters from user

            fprintf('\n=== Protection System Configuration ===\n\n');

            % Select relay type
            relayTypes = {'overcurrent', 'distance', 'differential', 'directional'};
            relayDescriptions = {
                'Overcurrent relay (50/51)', ...
                'Distance relay (21)', ...
                'Differential relay (87)', ...
                'Directional relay (67)'
            };

            fprintf('Available relay types:\n');
            relayIdx = prompt.getSelectOne(relayTypes, relayDescriptions, 'Select relay type:');
            obj.relayType = relayIdx;
            fprintf('Selected relay type: %s\n\n', obj.relayType);

            % System type
            systemTypes = {'distribution', 'transmission'};
            systemIdx = prompt.getMenuChoice(systemTypes, 'Select system type:');
            obj.systemType = systemTypes{systemIdx};

            % CT/VT ratios
            obj.CT_ratio = prompt.getNumericInputWithDefault( ...
                'Enter CT ratio (primary current): ', 400, [50, 10000]);

            obj.VT_ratio = prompt.getNumericInputWithDefault( ...
                'Enter VT ratio (primary voltage): ', 138000, [120, 500000]);

            % Pickup current (secondary)
            obj.pickupCurrent = prompt.getNumericInputWithDefault( ...
                'Enter pickup current (secondary, A): ', 5, [0.5, 20]);

            % Time dial setting
            obj.TDS = prompt.getNumericInputWithDefault( ...
                'Enter time dial setting (TDS): ', 1.0, [0.1, 15]);

            % Coordination time interval
            obj.CTI = prompt.getNumericInputWithDefault( ...
                'Enter coordination time interval (CTI, seconds): ', 0.3, [0.1, 0.5]);

            % Curve type
            curveTypes = {
                'IEEE moderately inverse', 'IEEE very inverse', 'IEEE extremely inverse', ...
                'IEC normal inverse', 'IEC very inverse', 'IEC extremely inverse'
            };
            curveIdx = prompt.getMenuChoice(curveTypes, 'Select time-current curve:');
            obj.curveType = curveTypes{curveIdx};
            fprintf('Selected curve: %s\n\n', obj.curveType);

            % Relay-specific parameters
            switch obj.relayType
                case 'distance'
                    obj.zones = obj.collectDistanceZones(prompt);
                case 'differential'
                    obj.slope = prompt.getNumericInputWithDefault( ...
                        'Enter differential slope (%%): ', 25, [10, 60]);
                case 'directional'
                    obj.maxTorqueAngle = prompt.getNumericInputWithDefault( ...
                        'Enter maximum torque angle (deg): ', 45, [30, 75]);
            end

            % Return parameters struct
            params.relayType = obj.relayType;
            params.CT_ratio = obj.CT_ratio;
            params.VT_ratio = obj.VT_ratio;
            params.pickupCurrent = obj.pickupCurrent;
            params.TDS = obj.TDS;
            params.CTI = obj.CTI;
            params.curveType = obj.curveType;
            params.systemType = obj.systemType;
            params.zones = obj.zones;
        end

        function zones = collectDistanceZones(obj, prompt)
            % Collect distance relay zone settings
            zones = struct('reach', {}, 'time', {}, 'angle', {});

            fprintf('\n--- Distance Relay Zone Settings ---\n');

            % Zone 1
            fprintf('Zone 1 (Instantaneous, 80-85%% of line):\n');
            zones(1).reach = prompt.getNumericInputWithDefault('  Reach (%% of line): ', 80, [50, 90]);
            zones(1).time = 0;  % Instantaneous
            zones(1).angle = prompt.getNumericInputWithDefault('  Characteristic angle (deg): ', 75, [60, 85]);

            % Zone 2
            fprintf('\nZone 2 (120-150%% with delay):\n');
            zones(2).reach = prompt.getNumericInputWithDefault('  Reach (%% of line): ', 120, [100, 200]);
            zones(2).time = prompt.getNumericInputWithDefault('  Time delay (s): ', 0.3, [0.2, 0.6]);
            zones(2).angle = zones(1).angle;

            % Zone 3
            fprintf('\nZone 3 (Backup protection):\n');
            zones(3).reach = prompt.getNumericInputWithDefault('  Reach (%% of line): ', 200, [150, 400]);
            zones(3).time = prompt.getNumericInputWithDefault('  Time delay (s): ', 1.0, [0.5, 2.0]);
            zones(3).angle = zones(1).angle;

            fprintf('\n');
        end

        function code = generateCode(obj, params)
            % Generate MATLAB code for protection analysis
            switch obj.relayType
                case 'overcurrent'
                    code = obj.generateOvercurrentCode();
                case 'distance'
                    code = obj.generateDistanceCode();
                case 'differential'
                    code = obj.generateDifferentialCode();
                case 'directional'
                    code = obj.generateDirectionalCode();
                otherwise
                    code = obj.generateOvercurrentCode();
            end
        end

        function script = generateOvercurrentCode(obj)
            % Generate overcurrent relay coordination code

            script = [
                '%% Overcurrent Relay Coordination' newline
                '% Auto-generated MATLAB code for relay coordination study' newline
                newline
                'clear; clc; close all;' newline
                newline
                '%% Relay Parameters' newline
                sprintf('CT_ratio = %d;' newline, obj.CT_ratio)
                sprintf('pickupSecondary = %g;  %% A (secondary)' newline, obj.pickupCurrent)
                sprintf('TDS = %g;  %% Time dial setting' newline, obj.TDS)
                sprintf('CTI = %g;  %% Coordination time interval (s)' newline, obj.CTI)
                sprintf('curveType = '''', obj.curveType)
                ''';' newline
                newline
                '%% Calculate Pickup Current (Primary)' newline
                'pickupPrimary = pickupSecondary * CT_ratio / 5;  %% A' newline
                'fprintf(''Pickup Current: %.1f A (primary)\n'', pickupPrimary);' newline
                newline
                '%% Time-Current Characteristic Curves' newline
                '% IEEE/IEC standard curves for overcurrent relays' newline
                newline
                '% Current multiple of pickup' newline
                'M = linspace(1, 40, 500);  %% Multiple of pickup current' newline
                newline
                '% IEEE curves coefficients' newline
                'switch curveType' newline
                '    case ''IEEE moderately inverse''' newline
                '        A = 0.17; B = 0.02; C = 0.018; D = 0.33; E = 5.95;' newline
                '        curveName = ''IEEE Moderately Inverse'';' newline
                '    case ''IEEE very inverse''' newline
                '        A = 19.61; B = 0.49; C = 0.04; D = 0.54; E = 22.06;' newline
                '        curveName = ''IEEE Very Inverse'';' newline
                '    case ''IEEE extremely inverse''' newline
                '        A = 28.2; B = 0.12; C = 0.03; D = 0.86; E = 36.0;' newline
                '        curveName = ''IEEE Extremely Inverse'';' newline
                '    case ''IEC normal inverse''' newline
                '        A = 0.14; B = 0; C = 0; D = 0.02; E = 1;' newline
                '        curveName = ''IEC Normal Inverse'';' newline
                '    case ''IEC very inverse''' newline
                '        A = 13.5; B = 0; C = 0; D = 0.25; E = 1;' newline
                '        curveName = ''IEC Very Inverse'';' newline
                '    case ''IEC extremely inverse''' newline
                '        A = 80; B = 0; C = 0; D = 0.14; E = 2;' newline
                '        curveName = ''IEC Extremely Inverse'';' newline
                '    otherwise' newline
                '        A = 0.17; B = 0.02; C = 0.018; D = 0.33; E = 5.95;' newline
                '        curveName = ''IEEE Moderately Inverse'';' newline
                'end' newline
                newline
                '%% Calculate Operating Time' newline
                '% Standard IEEE/IEC equation: t = TDS * (A/(M^B - C) + D + E/(M-F)^G)' newline
                '% Simplified: t = TDS * (A/(M^p - 1) + D)' newline
                newline
                't = zeros(size(M));' newline
                'for i = 1:length(M)' newline
                '    if M(i) > 1' newline
                sprintf('        t(i) = TDS * (%g./(M(i).^%g - 1) + %g);\n', ...
                    obj.getCurveCoefficients(obj.curveType))
                '    else' newline
                '        t(i) = inf;  %% Below pickup' newline
                '    end' newline
                'end' newline
                newline
            ];

            % Add coordination example
            script = [script
                '%% Relay Coordination Example' newline
                '% Calculate operating times for different fault locations' newline
                'faultCurrents = [500, 1000, 2000, 5000];  %% A (primary)' newline
                'fprintf(''\n--- Relay Operating Times ---\n'');' newline
                'fprintf(''Fault Current    Time Multiple    Operating Time\n'');' newline
                'fprintf(''------------------------------------------------\n'');' newline
                'for i = 1:length(faultCurrents)' newline
                '    If = faultCurrents(i);' newline
                '    M_val = If / pickupPrimary;' newline
                '    if M_val > 1' newline
                '        t_op = TDS * (A/(M_val^p - 1) + D);' newline
                '        fprintf(''%8d A         %8.2f          %8.3f s\n'', If, M_val, t_op);' newline
                '    else' newline
                '        fprintf(''%8d A         %8.2f          No trip\n'', If, M_val);' newline
                '    end' newline
                'end' newline
                newline
            ];

            % Add TCC plot
            script = [script
                '%% Time-Current Characteristic (TCC) Plot' newline
                'figure;' newline
                'loglog(M, t, ''b-'', ''LineWidth'', 2);' newline
                'hold on;' newline
                newline
                '% Plot for different TDS values' newline
                'TDS_values = [0.5, 1.0, 1.5, 2.0];' newline
                'colors = lines(length(TDS_values));' newline
                'for i = 1:length(TDS_values)' newline
                '    t_temp = zeros(size(M));' newline
                '    for j = 1:length(M)' newline
                '        if M(j) > 1' newline
                '            t_temp(j) = TDS_values(i) * (A/(M(j)^p - 1) + D);' newline
                '        else' newline
                '            t_temp(j) = inf;' newline
                '        end' newline
                '    end' newline
                '    loglog(M, t_temp, ''Color'', colors(i,:), ''LineWidth'', 1.5);' newline
                'end' newline
                newline
                'xlabel(''Current Multiple of Pickup'');' newline
                'ylabel(''Operating Time (seconds)'');' newline
                'title([''TCC Curve - '' curveName]);' newline
                'legend(arrayfun(@(x) sprintf(''TDS = %.1f'', x), TDS_values, ''UniformOutput'', false));' newline
                'grid on;' newline
                'ylim([0.01, 100]);' newline
                newline
                '%% Add Operating Regions' newline
                'xline(1, ''k--'', ''Pickup'', ''LineWidth'', 1);' newline
                'yline(CTI, ''r--'', ''CTI'', ''LineWidth'', 1);' newline
            ];

            % Add coordination calculation
            script = [script
                newline
                '%% Coordination Check' newline
                '% Compare primary and backup relay operating times' newline
                'TDS_primary = TDS;' newline
                'TDS_backup = TDS + CTI / (A/(2^p - 1) + D);  %% Simplified coordination' newline
                newline
                'fprintf(''\n--- Coordination Check ---\n'');' newline
                'fprintf(''Primary Relay TDS: %.2f\n'', TDS_primary);' newline
                'fprintf(''Backup Relay TDS: %.2f\n'', TDS_backup);' newline
                'fprintf(''Coordination Time Interval: %.3f s\n'', CTI);' newline
            ];
        end

        function script = generateDistanceCode(obj)
            % Generate distance relay protection code

            script = [
                '%% Distance Relay Protection' newline
                '% Auto-generated MATLAB code for distance relay analysis' newline
                newline
                'clear; clc; close all;' newline
                newline
                '%% System Parameters' newline
                sprintf('Z_line = 0.02 + 0.08i;  %% pu/mile (example)' newline)
                sprintf('lineLength = 50;  %% miles' newline)
                sprintf('CT_ratio = %d;' newline, obj.CT_ratio)
                sprintf('VT_ratio = %d;' newline, obj.VT_ratio)
                newline
                '%% Line Impedance' newline
                'Z_primary = Z_line * lineLength;  %% Total line impedance' newline
                'Z_secondary = Z_primary * (CT_ratio / VT_ratio);  %% Referred to secondary' newline
                newline
                '%% Distance Relay Zone Settings' newline
            ];

            % Add zone settings
            if isempty(obj.zones)
                script = [script
                    '% Zone 1: 80%% instantaneous' newline
                    'zone1_reach = 0.80;' newline
                    'zone1_time = 0;  %% Instantaneous' newline
                    newline
                    '% Zone 2: 120%% with 0.3s delay' newline
                    'zone2_reach = 1.20;' newline
                    'zone2_time = 0.3;  %% seconds' newline
                    newline
                    '% Zone 3: 200%% with 1.0s delay (backup)' newline
                    'zone3_reach = 2.00;' newline
                    'zone3_time = 1.0;  %% seconds' newline
                ];
            else
                for i = 1:length(obj.zones)
                    script = [script
                        sprintf('%% Zone %d Settings\n', i)
                        sprintf('zone%d_reach = %.2f;  %% %% of line\n', i, obj.zones(i).reach/100)
                        sprintf('zone%d_time = %.3f;  %% seconds\n', i, obj.zones(i).time)
                        newline
                    ];
                end
            end

            script = [script
                newline
                '%% Zone Impedance Settings (Secondary Ohms)' newline
                'Z1_zone = zone1_reach * Z_secondary;' newline
                'Z2_zone = zone2_reach * Z_secondary;' newline
                'Z3_zone = zone3_reach * Z_secondary;' newline
                newline
                'fprintf(''\n--- Distance Relay Zone Settings ---\n'');' newline
                'fprintf(''Line Impedance: %.4f /_ %.2f deg ohms (primary)\n'', ...' newline
                '    abs(Z_primary), angle(Z_primary)*180/pi);' newline
                'fprintf(''\nZone 1: %.4f ohms, %.0f%% reach, %s\n'', ...' newline
                '    abs(Z1_zone), zone1_reach*100, ...' newline
                '    ternary(zone1_time==0, ''instantaneous'', sprintf(''%.3f s'', zone1_time)));' newline
                'fprintf(''Zone 2: %.4f ohms, %.0f%% reach, %.3f s\n'', ...' newline
                '    abs(Z2_zone), zone2_reach*100, zone2_time);' newline
                'fprintf(''Zone 3: %.4f ohms, %.0f%% reach, %.3f s\n'', ...' newline
                '    abs(Z3_zone), zone3_reach*100, zone3_time);' newline
                newline
                '%% Impedance Diagram (R-X Plane)' newline
                'figure;' newline
                'theta = linspace(0, pi/2, 100);' newline
                newline
                '% Plot zone circles' newline
                'plot(abs(Z1_zone)*cos(theta), abs(Z1_zone)*sin(theta), ''r-'', ''LineWidth'', 2);' newline
                'hold on;' newline
                'plot(abs(Z2_zone)*cos(theta), abs(Z2_zone)*sin(theta), ''b-'', ''LineWidth'', 2);' newline
                'plot(abs(Z3_zone)*cos(theta), abs(Z3_zone)*sin(theta), ''g-'', ''LineWidth'', 2);' newline
                newline
                '% Plot line impedance' newline
                'plot([0, real(Z_secondary)], [0, imag(Z_secondary)], ''k-'', ''LineWidth'', 2);' newline
                'plot(real(Z_secondary), imag(Z_secondary), ''ko'', ''MarkerFaceColor'', ''k'');' newline
                newline
                'xlabel(''Resistance (ohms)'');' newline
                'ylabel(''Reactance (ohms)'');' newline
                'title(''Distance Relay Zone Characteristics'');' newline
                'legend(''Zone 1'', ''Zone 2'', ''Zone 3'', ''Line Impedance'');' newline
                'grid on;' newline
                'axis equal;' newline
                newline
                '%% Fault Location Calculation' newline
                '% Calculate apparent impedance for faults at different locations' newline
                'faultLocations = linspace(0, 1.5, 100);  %% 0 to 150%% of line' newline
                'Z_fault = faultLocations * Z_secondary;' newline
                newline
                '% Determine which zone would operate' newline
                'zone1_fault = abs(Z_fault) <= abs(Z1_zone);' newline
                'zone2_fault = abs(Z_fault) <= abs(Z2_zone) & ~zone1_fault;' newline
                'zone3_fault = abs(Z_fault) <= abs(Z3_zone) & ~zone2_fault & ~zone1_fault;' newline
                newline
                'figure;' newline
                'plot(faultLocations*100, abs(Z_fault), ''b-'', ''LineWidth'', 2);' newline
                'hold on;' newline
                'yline(abs(Z1_zone), ''r--'', ''Zone 1'');' newline
                'yline(abs(Z2_zone), ''b--'', ''Zone 2'');' newline
                'yline(abs(Z3_zone), ''g--'', ''Zone 3'');' newline
                'xlabel(''Fault Location (%% of line length)'');' newline
                'ylabel(''Apparent Impedance (ohms)'');' newline
                'title(''Fault Location vs Apparent Impedance'');' newline
                'grid on;' newline
            ];
        end

        function script = generateDifferentialCode(obj)
            % Generate differential relay protection code

            script = [
                '%% Differential Relay Protection' newline
                '% Auto-generated MATLAB code for differential protection' newline
                newline
                'clear; clc; close all;' newline
                newline
                '%% Differential Relay Parameters' newline
                sprintf('slope = %g;  %% Percentage slope' newline, obj.slope)
                sprintf('I_pickup = %g;  %% Minimum pickup current (A)' newline, obj.pickupCurrent)
                sprintf('CT_ratio = %d;' newline, obj.CT_ratio)
                newline
                '%% Differential Relay Operating Principle' newline
                '% I_diff = |I1 + I2| (differential current)' newline
                '% I_rest = |I1| + |I2| (restraint current)' newline
                '% Trip when: I_diff > slope * I_rest AND I_diff > I_pickup' newline
                newline
                '%% Generate Operating Characteristic' newline
                'I_rest = linspace(0, 10, 500);  %% Restraint current' newline
                'I_diff_trip = max(I_pickup, slope/100 * I_rest);  %% Trip boundary' newline
                newline
                'figure;' newline
                'plot(I_rest, I_diff_trip, ''r-'', ''LineWidth'', 2);' newline
                'fill([0, I_rest(end), I_rest(end), 0], ...' newline
                '    [I_pickup, I_diff_trip(2:end), 0, 0], ''r'', ''FaceAlpha'', 0.2);' newline
                'hold on;' newline
                newline
                '% Plot operating and restraining regions' newline
                'plot(I_rest, zeros(size(I_rest)), ''b--'');' newline
                'yline(I_pickup, ''k--'', ''Min Pickup'');' newline
                newline
                'xlabel(''Restraint Current (pu)'');' newline
                'ylabel(''Differential Current (pu)'');' newline
                'title([''Differential Relay Characteristic (Slope = '' num2str(slope) ''%%)'']);' newline
                'legend(''Trip Boundary'', ''Operate Region'', ''Min Pickup'');' newline
                'grid on;' newline
                newline
                '%% Example: Internal Fault' newline
                'I1_internal = 5;  %% Primary current' newline
                'I2_internal = -4;  %% Secondary current (reversed for internal fault)' newline
                'I_diff_internal = abs(I1_internal + I2_internal);' newline
                'I_rest_internal = abs(I1_internal) + abs(I2_internal);' newline
                'I_trip_boundary = max(I_pickup, slope/100 * I_rest_internal);' newline
                newline
                'fprintf(''\n--- Internal Fault Example ---\n'');' newline
                'fprintf(''I1 = %.2f A, I2 = %.2f A\n'', I1_internal, I2_internal);' newline
                'fprintf(''I_diff = %.2f A, I_rest = %.2f A\n'', I_diff_internal, I_rest_internal);' newline
                'fprintf(''Trip boundary = %.2f A\n'', I_trip_boundary);' newline
                'if I_diff_internal > I_trip_boundary' newline
                '    fprintf(''Result: TRIP\n'');' newline
                'else' newline
                '    fprintf(''Result: No trip\n'');' newline
                'end' newline
                newline
                '%% Example: External Fault (CT Saturation)' newline
                'I1_external = 10;  %% Primary current' newline
                'I2_external = -9;  %% Secondary (CT saturation causes mismatch)' newline
                'I_diff_external = abs(I1_external + I2_external);' newline
                'I_rest_external = abs(I1_external) + abs(I2_external);' newline
                'I_trip_boundary_ext = max(I_pickup, slope/100 * I_rest_external);' newline
                newline
                'fprintf(''\n--- External Fault (CT Saturation) ---\n'');' newline
                'fprintf(''I1 = %.2f A, I2 = %.2f A\n'', I1_external, I2_external);' newline
                'fprintf(''I_diff = %.2f A, I_rest = %.2f A\n'', I_diff_external, I_rest_external);' newline
                'fprintf(''Trip boundary = %.2f A\n'', I_trip_boundary_ext);' newline
                'if I_diff_external > I_trip_boundary_ext' newline
                '    fprintf(''Result: TRIP (misoperation!)\n'');' newline
                'else' newline
                '    fprintf(''Result: No trip (correct)\n'');' newline
                'end' newline
            ];
        end

        function script = generateDirectionalCode(obj)
            % Generate directional relay protection code

            script = [
                '%% Directional Overcurrent Relay' newline
                '% Auto-generated MATLAB code for directional protection' newline
                newline
                'clear; clc; close all;' newline
                newline
                '%% Directional Relay Parameters' newline
                sprintf('MTA = %d;  %% Maximum Torque Angle (degrees)' newline, obj.maxTorqueAngle)
                sprintf('pickupCurrent = %g;  %% A' newline, obj.pickupCurrent)
                newline
                '%% Directional Element Operating Principle' newline
                '% Operates when current flows in forward direction' newline
                '% Characteristic angle defines the reference' newline
                newline
                '%% Generate Directional Characteristic' newline
                'theta = linspace(0, 2*pi, 360);' newline
                'MTA_rad = MTA * pi / 180;' newline
                newline
                '% Operating region: within +/- 90 degrees of MTA' newline
                'forward_region = mod(theta - MTA_rad + pi, 2*pi) - pi;' newline
                'operating = abs(forward_region) < pi/2;' newline
                newline
                'figure;' newline
                'polarplot(theta, ones(size(theta)), ''k--'');' newline
                'hold on;' newline
                newline
                '% Plot operating region' newline
                'theta_op = theta(operating);' newline
                'polarplot(theta_op, 0.5 * ones(size(theta_op)), ''r'', ''LineWidth'', 3);' newline
                newline
                '% Plot MTA line' newline
                'polarplot([MTA_rad, MTA_rad + pi], [0, 1], ''b-'', ''LineWidth'', 2);' newline
                newline
                'title([''Directional Characteristic (MTA = '' num2str(MTA) ''°)'']);' newline
                newline
                '%% Phase Angle Comparison' newline
                '% Example: Calculate operating torque for different fault angles' newline
                'faultAngles = -180:10:180;  %% degrees' newline
                'torque = cosd(faultAngles - MTA);' newline
                newline
                'figure;' newline
                'plot(faultAngles, torque, ''b-'', ''LineWidth'', 2);' newline
                'hold on;' newline
                'yline(0, ''k--'');' newline
                'xline(MTA, ''r--'', ''MTA'');' newline
                'xlabel(''Fault Current Angle (degrees)'');' newline
                'ylabel(''Relative Torque'');' newline
                'title(''Directional Element Torque vs Fault Angle'');' newline
                'grid on;' newline
                newline
                'fprintf(''\n--- Directional Relay Operation ---\n'');' newline
                'for i = 1:length(faultAngles)' newline
                '    if torque(i) > 0' newline
                '        direction = ''Forward'';' newline
                '    else' newline
                '        direction = ''Reverse'';' newline
                '    end' newline
                '    if mod(i,10) == 1  %% Print every 10th value' newline
                '        fprintf(''Angle %4d°: Torque = %6.3f, %s\n'', ...' newline
                '            faultAngles(i), torque(i), direction);' newline
                '    end' newline
                'end' newline
            ];
        end

        function [A, p, D] = getCurveCoefficients(obj, curveType)
            % Get IEEE curve coefficients
            switch curveType
                case 'IEEE moderately inverse'
                    A = 0.17; p = 0.02; D = 0.33;
                case 'IEEE very inverse'
                    A = 19.61; p = 0.49; D = 0.54;
                case 'IEEE extremely inverse'
                    A = 28.2; p = 0.12; D = 0.86;
                otherwise
                    A = 0.17; p = 0.02; D = 0.33;
            end
        end

        function generateSimulinkModel(obj, params)
            % Generate Simulink model for relay simulation
            fprintf('\nGenerating Simulink relay model...\n\n');

            script = [
                '% Relay Simulation Model' newline
                'modelName = ''relay_simulation'';' newline
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
                '% Add fault block' newline
                'add_block(''powerlib/Faults/Three-Phase Fault'', ...' newline
                '    [modelName ''/Fault'']);' newline
                newline
                '% Add relay model' newline
                'add_block(''simulink/Logic and Bit Operations/Relay'', ...' newline
                '    [modelName ''/Overcurrent Relay'']);' newline
                newline
                '% Add trip logic' newline
                'add_block(''simulink/Sinks/Scope'', [modelName ''/Scope'']);' newline
                newline
                'open_system(modelName);' newline
                'fprintf(''Relay simulation model created: %s.slx\n'', modelName);' newline
            ];

            fprintf('%s\n', script);
        end

        function handleInput(obj, userInput)
            % Handle contextual user input
            fprintf('Protection module received: %s\n', userInput);
        end

        function str = ternary(condition, trueVal, falseVal)
            % Helper function for ternary operation
            if condition
                str = trueVal;
            else
                str = falseVal;
            end
        end
    end
end
