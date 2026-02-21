classdef Validator < handle
    % Validator - Validates parameters against engineering constraints
    %
    % Ensures that user-provided parameters are within acceptable ranges
    % and consistent with power system engineering practices.

    methods
        function obj = Validator()
            % Constructor
        end

        function result = validateParameters(obj, params, domain)
            % Validate parameters for a specific domain
            %
            % Returns struct with isValid flag and warnings array

            result.isValid = true;
            result.warnings = {};
            result.errors = {};

            switch domain
                case 'power_flow'
                    result = obj.validatePowerFlowParams(params, result);

                case 'fault_analysis'
                    result = obj.validateFaultParams(params, result);

                case 'protection'
                    result = obj.validateProtectionParams(params, result);

                case 'stability'
                    result = obj.validateStabilityParams(params, result);

                case 'renewable_integration'
                    result = obj.validateRenewableParams(params, result);
            end

            % Check for any errors
            if ~isempty(result.errors)
                result.isValid = false;
            end
        end

        function result = validatePowerFlowParams(obj, params, result)
            % Validate power flow parameters

            % Check number of buses
            if isfield(params, 'numBuses')
                if params.numBuses < 2
                    result.errors{end+1} = 'At least 2 buses required for power flow';
                end
                if params.numBuses > 10000
                    result.warnings{end+1} = 'Large system: Consider sparse matrix techniques';
                end
            end

            % Check tolerance
            if isfield(params, 'tolerance')
                if params.tolerance <= 0
                    result.errors{end+1} = 'Convergence tolerance must be positive';
                elseif params.tolerance > 1e-3
                    result.warnings{end+1} = 'Large tolerance may give inaccurate results';
                elseif params.tolerance < 1e-12
                    result.warnings{end+1} = 'Very small tolerance may cause convergence issues';
                end
            end

            % Check max iterations
            if isfield(params, 'maxIter')
                if params.maxIter < 1
                    result.errors{end+1} = 'Maximum iterations must be at least 1';
                elseif params.maxIter > 1000
                    result.warnings{end+1} = 'Large iteration limit may indicate convergence problems';
                end
            end

            % Check method
            if isfield(params, 'method')
                validMethods = {'Newton-Raphson', 'Gauss-Seidel', 'Fast Decoupled', 'DC Power Flow'};
                if ~any(strcmp(params.method, validMethods))
                    result.errors{end+1} = sprintf('Invalid method. Choose from: %s', strjoin(validMethods, ', '));
                end
            end

            % Check for slack bus
            if isfield(params, 'busData')
                hasSlack = false;
                for i = 1:length(params.busData)
                    if isfield(params.busData(i), 'type') && params.busData(i).type == 3
                        hasSlack = true;
                        break;
                    end
                end
                if ~hasSlack
                    result.errors{end+1} = 'Power flow requires at least one slack bus';
                end
            end
        end

        function result = validateFaultParams(obj, params, result)
            % Validate fault analysis parameters

            % Check fault impedance
            if isfield(params, 'faultImpedance')
                if params.faultImpedance < 0
                    result.errors{end+1} = 'Fault impedance cannot be negative';
                elseif params.faultImpedance > 1
                    result.warnings{end+1} = 'Very high fault impedance may indicate unrealistic scenario';
                end
            end

            % Check X/R ratio
            if isfield(params, 'XR_ratio')
                if params.XR_ratio < 0.1
                    result.warnings{end+1} = 'Very low X/R ratio - typical range is 1-15';
                elseif params.XR_ratio > 50
                    result.warnings{end+1} = 'Very high X/R ratio - typical range is 1-15';
                end
            end

            % Check fault type
            if isfield(params, 'faultType')
                validTypes = {'three_phase', 'single_line_ground', 'line_to_line', 'double_line_ground'};
                if ~any(strcmp(params.faultType, validTypes))
                    result.errors{end+1} = sprintf('Invalid fault type. Choose from: %s', strjoin(validTypes, ', '));
                end
            end

            % Check voltage level
            if isfield(params, 'voltageLevel')
                if params.voltageLevel < 0.12
                    result.warnings{end+1} = 'Very low voltage - verify base kV';
                elseif params.voltageLevel > 1000
                    result.warnings{end+1} = 'EHV/UHV system - consider special phenomena';
                end
            end
        end

        function result = validateProtectionParams(obj, params, result)
            % Validate protection system parameters

            % Check CTI (Coordination Time Interval)
            if isfield(params, 'CTI')
                if params.CTI < 0.1
                    result.warnings{end+1} = 'CTI < 0.1s may cause coordination issues';
                elseif params.CTI > 0.5
                    result.warnings{end+1} = 'Large CTI increases fault clearing time';
                end
            end

            % Check TDS (Time Dial Setting)
            if isfield(params, 'TDS')
                if params.TDS < 0.1
                    result.warnings{end+1} = 'Very low TDS - check relay limits';
                elseif params.TDS > 15
                    result.warnings{end+1} = 'Very high TDS - may cause slow operation';
                end
            end

            % Check CT ratio
            if isfield(params, 'CT_ratio')
                if params.CT_ratio <= 1
                    result.errors{end+1} = 'Invalid CT ratio';
                end
            end

            % Check pickup current
            if isfield(params, 'pickupCurrent')
                if params.pickupCurrent <= 0
                    result.errors{end+1} = 'Pickup current must be positive';
                end
            end

            % Check relay type
            if isfield(params, 'relayType')
                validTypes = {'overcurrent', 'distance', 'differential', 'directional'};
                if ~any(strcmp(params.relayType, validTypes))
                    result.errors{end+1} = sprintf('Invalid relay type. Choose from: %s', strjoin(validTypes, ', '));
                end
            end
        end

        function result = validateStabilityParams(obj, params, result)
            % Validate stability analysis parameters

            % Check simulation time step
            if isfield(params, 'stepSize')
                if params.stepSize <= 0
                    result.errors{end+1} = 'Step size must be positive';
                elseif params.stepSize > 0.1
                    result.warnings{end+1} = 'Large step size may miss fast dynamics';
                elseif params.stepSize < 1e-6
                    result.warnings{end+1} = 'Very small step size increases computation time';
                end
            end

            % Check simulation time
            if isfield(params, 'simulationTime')
                if params.simulationTime <= 0
                    result.errors{end+1} = 'Simulation time must be positive';
                elseif params.simulationTime < 1
                    result.warnings{end+1} = 'Short simulation may not capture slow dynamics';
                elseif params.simulationTime > 100
                    result.warnings{end+1} = 'Long simulation time - consider output interval';
                end
            end

            % Check study type
            if isfield(params, 'studyType')
                validTypes = {'transient', 'small_signal', 'voltage', 'frequency'};
                if ~any(strcmp(params.studyType, validTypes))
                    result.errors{end+1} = sprintf('Invalid study type. Choose from: %s', strjoin(validTypes, ', '));
                end
            end

            % Check solver
            if isfield(params, 'solver')
                validSolvers = {'ode45', 'ode23t', 'ode15s', 'ode23'};
                if ~any(strcmp(params.solver, validSolvers))
                    result.warnings{end+1} = sprintf('Consider using ode23t for power system stability');
                end
            end
        end

        function result = validateRenewableParams(obj, params, result)
            % Validate renewable integration parameters

            % Check capacity
            if isfield(params, 'capacity')
                if params.capacity <= 0
                    result.errors{end+1} = 'Capacity must be positive';
                end
            end

            % Check penetration level
            if isfield(params, 'penetrationLevel')
                if params.penetrationLevel < 0
                    result.errors{end+1} = 'Penetration level cannot be negative';
                elseif params.penetrationLevel > 100
                    result.errors{end+1} = 'Penetration level cannot exceed 100%';
                elseif params.penetrationLevel > 50
                    result.warnings{end+1} = 'High penetration - consider system strength and stability';
                end
            end

            % Check resource type
            if isfield(params, 'resourceType')
                validTypes = {'solar_pv', 'wind', 'storage', 'hybrid'};
                if ~any(strcmp(params.resourceType, validTypes))
                    result.errors{end+1} = sprintf('Invalid resource type. Choose from: %s', strjoin(validTypes, ', '));
                end
            end

            % Check control mode
            if isfield(params, 'controlMode')
                validModes = {'PQ', 'PV', 'droop', 'grid_forming'};
                if ~any(strcmp(params.controlMode, validModes))
                    result.warnings{end+1} = 'Verify control mode is appropriate for study type';
                end
            end
        end

        function isValid = checkVoltageRange(obj, voltage, baseVoltage)
            % Check if voltage is within acceptable range
            tolerance = 0.05;  % 5% tolerance
            minVoltage = baseVoltage * (1 - tolerance);
            maxVoltage = baseVoltage * (1 + tolerance);

            isValid = (voltage >= minVoltage) && (voltage <= maxVoltage);
        end

        function isValid = checkThermalLimit(obj, current, rating)
            % Check if current is within thermal rating
            safetyFactor = 0.9;  % 90% of rating
            isValid = current <= (rating * safetyFactor);
        end

        function isValid = checkStabilityLimit(obj, powerTransfer, maxTransfer)
            % Check if power transfer is within stability limit
            margin = 0.8;  % 80% of max for safety margin
            isValid = powerTransfer <= (maxTransfer * margin);
        end
    end
end
