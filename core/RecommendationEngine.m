classdef RecommendationEngine < handle
    % RecommendationEngine - Provides best practice recommendations
    %
    % This engine contains domain-specific knowledge about power system
    % engineering best practices, standard parameters, and IEEE/IEC guidelines.

    properties
        % Standard parameter ranges from IEEE/IEC standards
        transmissionXR       = [10, 15];      % Transmission X/R ratio
        distributionXR       = [1, 5];        % Distribution X/R ratio
        coordinationTimeInt  = [0.2, 0.4];    % CTI in seconds
        faultImpedance       = [0.01, 0.1];   % Fault impedance (pu)
        stabilityStep        = [0.001, 0.01]; % Simulation step for stability
        powerFlowTolerance   = 1e-6;          % Convergence tolerance
        voltageTolerance     = [0.95, 1.05];  % Normal voltage range (pu)
        frequencyTolerance   = [59.5, 60.5];  % Normal frequency range (Hz)

        % Solver recommendations
        solvers              = {'ode23t', 'ode15s', 'ode45'};
        recommendedSolver    = 'ode23t';      % Default for power systems
    end

    methods
        function obj = RecommendationEngine()
            % Constructor - initialize recommendation database
        end

        function recommendations = getRecommendations(obj, domain, parameters)
            % Get domain-specific recommendations
            recommendations = {};

            switch domain
                case 'power_flow'
                    recommendations = obj.getPowerFlowRecommendations(parameters);

                case 'fault_analysis'
                    recommendations = obj.getFaultRecommendations(parameters);

                case 'protection'
                    recommendations = obj.getProtectionRecommendations(parameters);

                case 'stability'
                    recommendations = obj.getStabilityRecommendations(parameters);

                case 'renewable_integration'
                    recommendations = obj.getRenewableRecommendations(parameters);

                otherwise
                    recommendations = obj.getGeneralRecommendations();
            end
        end

        function recs = getPowerFlowRecommendations(obj, params)
            % Power flow specific recommendations
            recs = {
                'Use Newton-Raphson method for large systems (>50 buses) for faster convergence',
                'Gauss-Seidel is suitable for small systems and educational purposes',
                sprintf('Set convergence tolerance to %.0e for accurate results', obj.powerFlowTolerance),
                'Ensure slack bus has sufficient capacity to balance the system',
                'Check voltage magnitudes are within acceptable range (0.95-1.05 pu)',
                'For ill-conditioned systems, try Fast Decoupled method with XB version'
            };

            % Add specific recommendations based on parameters
            if isfield(params, 'numBuses') && params.numBuses > 100
                recs{end+1} = 'Consider using sparse matrix techniques for large systems';
            end

            if isfield(params, 'method') && strcmp(params.method, 'Gauss-Seidel')
                if isfield(params, 'numBuses') && params.numBuses > 50
                    recs{end+1} = 'Warning: Gauss-Seidel may converge slowly for large systems. Consider Newton-Raphson.';
                end
            end
        end

        function recs = getFaultRecommendations(obj, params)
            % Fault analysis specific recommendations
            recs = {
                sprintf('Typical fault impedance range: %.2f to %.2f pu', obj.faultImpedance(1), obj.faultImpedance(2)),
                'For transmission systems, use X/R ratio between 10-15',
                'For distribution systems, use X/R ratio between 1-5',
                'Three-phase faults typically produce the highest fault currents',
                'Single line-to-ground faults are most common in practice',
                'Consider DC offset in breaker rating calculations'
            };

            % Add fault-specific recommendations
            if isfield(params, 'faultType')
                switch params.faultType
                    case 'three_phase'
                        recs{end+1} = 'Three-phase fault: Use for breaker rating and worst-case analysis';
                    case 'single_line_ground'
                        recs{end+1} = 'SLG fault: Most common (70-80% of faults). Consider ground resistance.';
                    case 'line_to_line'
                        recs{end+1} = 'LL fault: Typically occurs during fault clearing or lightning strikes';
                    case 'double_line_ground'
                        recs{end+1} = 'DLG fault: More severe than LL but less than three-phase';
                end
            end

            if isfield(params, 'voltageLevel') && params.voltageLevel >= 230
                recs{end+1} = 'High voltage system: Consider transient recovery voltage (TRV) studies';
            end
        end

        function recs = getProtectionRecommendations(obj, params)
            % Protection system specific recommendations
            recs = {
                sprintf('Coordination Time Interval (CTI): %.1f to %.1f seconds', ...
                    obj.coordinationTimeInt(1), obj.coordinationTimeInt(2)),
                'Typical CT saturation factor: 1.5 to 2.0 for relaying applications',
                'Use class C or K CTs for differential protection',
                'Distance relay Zone 1: Set to 80-85% of line length',
                'Distance relay Zone 2: Set to 120-150% of line length with 0.3-0.5s delay',
                'Overcurrent pickup: 125-150% of maximum load current'
            };

            % Add relay-specific recommendations
            if isfield(params, 'relayType')
                switch params.relayType
                    case 'overcurrent'
                        recs{end+1} = 'Time dial setting (TDS): Start with 0.5-1.5 for coordination';
                        recs{end+1} = 'Use IEC or IEEE inverse curves based on utility standard';
                    case 'distance'
                        recs{end+1} = 'Ensure proper coordination with adjacent line protections';
                        recs{end+1} = 'Consider power swing blocking for long lines';
                    case 'differential'
                        recs{end+1} = 'Set slope to 20-40% to account for CT errors and tap changing';
                        recs{end+1} = 'Use harmonic restraint to prevent operation on transformer inrush';
                end
            end

            if isfield(params, 'systemType') && strcmp(params.systemType, 'distribution')
                recs{end+1} = 'For distribution: Consider fuse coordination and recloser settings';
            end
        end

        function recs = getStabilityRecommendations(obj, params)
            % Stability analysis specific recommendations
            recs = {
                sprintf('Simulation time step: %.3f to %.3f seconds for dynamic studies', ...
                    obj.stabilityStep(1), obj.stabilityStep(2)),
                'Use trapezoidal integration (ode23t) for power system stability',
                'Simulation time: 10-20 seconds for transient stability assessment',
                'Include governor and AVR models for accurate long-term dynamics',
                'Consider load modeling: ZIP or induction motor models for accuracy',
                'Small signal analysis: Check eigenvalues for damping ratios > 5%'
            };

            % Add stability-specific recommendations
            if isfield(params, 'studyType')
                switch params.studyType
                    case 'transient'
                        recs{end+1} = 'Critical clearing time: Test multiple fault locations and clearing times';
                        recs{end+1} = 'Monitor rotor angles: Separation > 180 degrees indicates instability';
                    case 'small_signal'
                        recs{end+1} = 'Identify poorly damped modes (damping < 5%) for PSS tuning';
                        recs{end+1} = 'Consider participation factors to identify generator contributions';
                    case 'voltage'
                        recs{end+1} = 'P-V curves: Identify nose point for voltage collapse margin';
                        recs{end+1} = 'Consider reactive power limits and load tap changers';
                end
            end

            if isfield(params, 'simulationTime') && params.simulationTime < 5
                recs{end+1} = 'Warning: Short simulation time may not capture slow dynamics (governors, LTCs)';
            end
        end

        function recs = getRenewableRecommendations(obj, params)
            % Renewable integration specific recommendations
            recs = {
                'PV inverter: Use IEEE 1547 compliant models for grid interconnection',
                'Wind turbine: Type 3 (DFIG) or Type 4 (full converter) for modern turbines',
                'Include low voltage ride-through (LVRT) capability in inverter models',
                'Consider flicker and harmonic analysis for inverter-based resources',
                'Energy storage: Model SOC dynamics for realistic dispatch behavior',
                'Grid-forming inverters: Essential for microgrid islanded operation'
            };

            % Add renewable-specific recommendations
            if isfield(params, 'resourceType')
                switch params.resourceType
                    case 'solar_pv'
                        recs{end+1} = 'PV model: Include irradiance and temperature dependencies';
                        recs{end+1} = 'Inverter control: dq-frame or phasor model based on study type';
                    case 'wind'
                        recs{end+1} = 'Wind model: Include wind speed variations and turbine dynamics';
                        recs{end+1} = 'DFIG: Model rotor-side and grid-side converter controls';
                    case 'storage'
                        recs{end+1} = 'Battery model: Include internal resistance and capacity limits';
                        recs{end+1} = 'BMS: Model state of charge (SOC) and power limits';
                end
            end

            if isfield(params, 'penetrationLevel') && params.penetrationLevel > 30
                recs{end+1} = 'High penetration (>30%): Consider system strength and short circuit ratio';
                recs{end+1} = 'Verify frequency response and inertia requirements are met';
            end
        end

        function recs = getGeneralRecommendations(obj)
            % General power system engineering recommendations
            recs = {
                'Always verify model results against hand calculations for simple cases',
                'Use per-unit system for consistency across different voltage levels',
                'Document all assumptions and base values used in the analysis',
                'Validate model with field measurements when available',
                'Consider worst-case scenarios for design and protection coordination'
            };
        end

        function suggestion = suggestParameter(obj, paramName, context)
            % Suggest a standard value for a parameter
            suggestion = [];

            switch paramName
                case 'XR_ratio'
                    if isfield(context, 'systemType')
                        if strcmp(context.systemType, 'transmission')
                            suggestion = mean(obj.transmissionXR);
                        else
                            suggestion = mean(obj.distributionXR);
                        end
                    else
                        suggestion = 10;  % Default
                    end

                case 'CTI'
                    suggestion = mean(obj.coordinationTimeInt);

                case 'faultZ'
                    suggestion = mean(obj.faultImpedance);

                case 'stepSize'
                    suggestion = 0.001;

                case 'tolerance'
                    suggestion = obj.powerFlowTolerance;

                case 'voltageMin'
                    suggestion = obj.voltageTolerance(1);

                case 'voltageMax'
                    suggestion = obj.voltageTolerance(2);

                case 'solver'
                    suggestion = obj.recommendedSolver;
            end
        end

        function info = getStandardInfo(obj, standardName)
            % Get information about power system standards
            standards = struct();
            standards.IEEE_1547 = 'Interconnection and Interoperability of Distributed Energy Resources';
            standards.IEEE_399 = 'Brown Book - Power Systems Analysis (Industrial and Commercial)';
            standards.IEEE_241 = 'Gray Book - Electric Power Systems Design';
            standards.IEEE_141 = 'Red Book - Electric Power Distribution for Industrial Plants';
            standards.IEC_60255 = 'Measuring Relays and Protection Equipment';
            standards.IEC_61850 = 'Communication Networks and Systems for Power Utility Automation';

            if isfield(standards, standardName)
                info = standards.(standardName);
            else
                info = 'Standard not found in database';
            end
        end
    end
end
