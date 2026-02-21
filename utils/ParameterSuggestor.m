classdef ParameterSuggestor < handle
    % ParameterSuggestor - Provides smart default values for parameters
    %
    % Suggests appropriate default values based on system type,
    % voltage level, and engineering best practices.

    properties
        % Standard voltage levels (kV)
        standardVoltages = [0.12, 0.208, 0.24, 0.48, 2.4, 4.16, 13.8, 34.5, 69, 115, 138, 161, 230, 345, 500, 765];

        % Typical parameters by voltage level
        typicalXR = struct();
        typicalFaultZ = struct();
        typicalSCRatio = struct();
    end

    methods
        function obj = ParameterSuggestor()
            % Constructor - initialize typical parameters
            obj.typicalXR.distribution = 3;
            obj.typicalXR.subtransmission = 8;
            obj.typicalXR.transmission = 12;
            obj.typicalXR.EHV = 15;

            obj.typicalFaultZ.distribution = 0.1;
            obj.typicalFaultZ.subtransmission = 0.05;
            obj.typicalFaultZ.transmission = 0.02;
            obj.typicalFaultZ.EHV = 0.01;
        end

        function suggestion = suggestDefaults(obj, domain, context)
            % Suggest default parameters for a domain
            suggestion = struct();

            switch domain
                case 'power_flow'
                    suggestion = obj.suggestPowerFlowDefaults(context);

                case 'fault_analysis'
                    suggestion = obj.suggestFaultDefaults(context);

                case 'protection'
                    suggestion = obj.suggestProtectionDefaults(context);

                case 'stability'
                    suggestion = obj.suggestStabilityDefaults(context);

                case 'renewable_integration'
                    suggestion = obj.suggestRenewableDefaults(context);
            end
        end

        function defaults = suggestPowerFlowDefaults(obj, context)
            % Suggest defaults for power flow analysis
            defaults.method = 'Newton-Raphson';
            defaults.tolerance = 1e-6;
            defaults.maxIter = 50;
            defaults.baseMVA = 100;
            defaults.accelerationFactor = 1.6;  % For Gauss-Seidel

            % Voltage level specific
            if isfield(context, 'voltageLevel')
                if context.voltageLevel >= 230
                    defaults.XR_ratio = 12;
                elseif context.voltageLevel >= 69
                    defaults.XR_ratio = 10;
                else
                    defaults.XR_ratio = 5;
                end
            else
                defaults.XR_ratio = 10;
            end
        end

        function defaults = suggestFaultDefaults(obj, context)
            % Suggest defaults for fault analysis
            defaults.faultType = 'three_phase';
            defaults.baseMVA = 100;

            % Determine system type from voltage
            if isfield(context, 'voltageLevel')
                vLevel = context.voltageLevel;
                if vLevel >= 230
                    sysType = 'EHV';
                elseif vLevel >= 115
                    sysType = 'transmission';
                elseif vLevel >= 34.5
                    sysType = 'subtransmission';
                else
                    sysType = 'distribution';
                end

                defaults.XR_ratio = obj.typicalXR.(sysType);
                defaults.faultImpedance = obj.typicalFaultZ.(sysType);
            else
                defaults.XR_ratio = 10;
                defaults.faultImpedance = 0.02;
            end

            % Pre-fault voltage
            defaults.prefaultVoltage = 1.0;

            % Include DC offset
            defaults.includeDCoffset = true;
        end

        function defaults = suggestProtectionDefaults(obj, context)
            % Suggest defaults for protection coordination
            defaults.CTI = 0.3;  % Coordination Time Interval
            defaults.TDS = 1.0;  % Time Dial Setting
            defaults.curveType = 'IEEE moderately inverse';

            % CT/VT ratios based on voltage
            if isfield(context, 'voltageLevel')
                vLevel = context.voltageLevel;

                % Suggest CT ratio based on typical currents
                if vLevel >= 230
                    defaults.CT_ratio = 2000;
                    defaults.VT_ratio = 2000;  % 230kV/sqrt(3) : 115V
                elseif vLevel >= 115
                    defaults.CT_ratio = 1200;
                    defaults.VT_ratio = 1000;
                elseif vLevel >= 34.5
                    defaults.CT_ratio = 600;
                    defaults.VT_ratio = 300;
                else
                    defaults.CT_ratio = 400;
                    defaults.VT_ratio = 120;
                end
            else
                defaults.CT_ratio = 400;
                defaults.VT_ratio = 120;
            end

            % Pickup current (125% of typical load)
            if isfield(context, 'loadCurrent')
                defaults.pickupCurrent = 1.25 * context.loadCurrent;
            else
                defaults.pickupCurrent = 100;  % Default amps
            end

            % Instantaneous setting
            defaults.instantaneousMultiplier = 1.25;  % Above max fault
        end

        function defaults = suggestStabilityDefaults(obj, context)
            % Suggest defaults for stability analysis
            defaults.studyType = 'transient';
            defaults.simulationTime = 10;  % seconds
            defaults.stepSize = 0.001;  % 1ms
            defaults.solver = 'ode23t';

            % Disturbance settings
            defaults.disturbanceType = 'three_phase_fault';
            defaults.faultClearingTime = 0.1;  % 6 cycles at 60Hz

            % Generator modeling
            defaults.generatorModel = 'classical';  % or 'detailed'
            defaults.includeGovernor = false;
            defaults.includeAVR = false;
            defaults.includePSS = false;

            % Load modeling
            defaults.loadModel = 'constant_impedance';

            % For detailed studies
            if isfield(context, 'studyDetail') && strcmp(context.studyDetail, 'detailed')
                defaults.generatorModel = 'detailed';
                defaults.includeGovernor = true;
                defaults.includeAVR = true;
                defaults.stepSize = 0.0001;  % 0.1ms for detailed
            end
        end

        function defaults = suggestRenewableDefaults(obj, context)
            % Suggest defaults for renewable integration
            defaults.controlMode = 'PQ';
            defaults.powerFactor = 1.0;

            % Resource-specific defaults
            if isfield(context, 'resourceType')
                switch context.resourceType
                    case 'solar_pv'
                        defaults.irradiance = 1000;  % W/m^2
                        defaults.temperature = 25;   % Celsius
                        defaults.inverterEfficiency = 0.97;

                    case 'wind'
                        defaults.windSpeed = 12;     % m/s
                        defaults.turbineType = 'Type4';
                        defaults.rotorDiameter = 100; % meters

                    case 'storage'
                        defaults.capacity_kWh = 1000;
                        defaults.maxPower_kW = 250;
                        defaults.initialSOC = 0.5;
                        defaults.efficiency = 0.95;
                end
            end

            % Grid connection
            defaults.includeLVRT = true;
            defaults.includeFilter = true;
        end

        function voltage = suggestVoltageLevel(obj, application)
            % Suggest standard voltage level for application
            switch lower(application)
                case {'residential', 'house'}
                    voltage = 0.12;  % 120V
                case {'commercial', 'building'}
                    voltage = 0.208;  % 208V
                case {'industrial', 'factory'}
                    voltage = 4.16;  % 4.16kV
                case {'distribution', 'feeder'}
                    voltage = 13.8;  % 13.8kV
                case {'subtransmission', 'substation'}
                    voltage = 69;  % 69kV
                case {'transmission', 'trans'}
                    voltage = 138;  % 138kV
                case {'high_voltage', 'HV'}
                    voltage = 230;  % 230kV
                case {'EHV', 'extra_high'}
                    voltage = 500;  % 500kV
                otherwise
                    voltage = 13.8;  % Default
            end
        end

        function ctRatio = suggestCTRatio(obj, primaryCurrent)
            % Suggest CT ratio based on primary current
            % Ensures secondary current is 5A at rated conditions

            standardRatios = [50, 100, 150, 200, 300, 400, 600, 800, 1000, ...
                              1200, 1500, 2000, 2500, 3000, 4000, 5000];

            % Find smallest standard ratio that gives <= 5A secondary
            for i = 1:length(standardRatios)
                ratio = standardRatios(i);
                secondaryCurrent = primaryCurrent * 5 / ratio;
                if secondaryCurrent <= 5
                    ctRatio = ratio;
                    return;
                end
            end

            % If none found, use largest
            ctRatio = standardRatios(end);
        end

        function baseMVA = suggestBaseMVA(obj, systemMVA)
            % Suggest appropriate base MVA
            standardBases = [1, 10, 25, 50, 100, 200, 500, 1000];

            % Find closest standard base
            [~, idx] = min(abs(standardBases - systemMVA));
            baseMVA = standardBases(idx);
        end

        function params = getSystemTypeParams(obj, systemType)
            % Get typical parameters for system type
            switch lower(systemType)
                case 'distribution'
                    params.voltageRange = [0.12, 34.5];
                    params.XR_ratio = 3;
                    params.R1_ratio = 0.8;  % R1/X1
                    params.R0_ratio = 3.0;  % R0/X0
                    params.X0_ratio = 2.5;  % X0/X1

                case 'subtransmission'
                    params.voltageRange = [34.5, 69];
                    params.XR_ratio = 8;
                    params.R1_ratio = 0.3;
                    params.R0_ratio = 2.0;
                    params.X0_ratio = 2.2;

                case 'transmission'
                    params.voltageRange = [69, 230];
                    params.XR_ratio = 12;
                    params.R1_ratio = 0.1;
                    params.R0_ratio = 1.5;
                    params.X0_ratio = 2.0;

                case {'EHV', 'extra_high'}
                    params.voltageRange = [230, 765];
                    params.XR_ratio = 15;
                    params.R1_ratio = 0.05;
                    params.R0_ratio = 1.2;
                    params.X0_ratio = 1.8;

                otherwise
                    params.XR_ratio = 10;
                    params.R1_ratio = 0.2;
            end
        end
    end
end
