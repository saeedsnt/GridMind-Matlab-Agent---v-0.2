classdef PowerSystemClassifier < handle
    % PowerSystemClassifier - Categorizes power system engineering tasks
    %
    % Provides domain-specific information, parameters, and module mappings
    % for different power system analysis categories.

    methods
        function obj = PowerSystemClassifier()
            % Constructor
        end

        function domainInfo = getDomainInfo(obj, domainName)
            % Get information about a specific power system domain
            %
            % Returns struct with domain details or empty if not found

            domains = obj.getAllDomains();

            if isfield(domains, domainName)
                domainInfo = domains.(domainName);
            else
                domainInfo = [];
            end
        end

        function domains = getAllDomains(obj)
            % Get all available power system domains with their details

            % Power Flow Analysis
            domains.power_flow.name = 'Power Flow Analysis';
            domains.power_flow.description = 'Load flow studies for steady-state analysis';
            domains.power_flow.module = 'PowerFlowModule';
            domains.power_flow.implementation = 'code';
            domains.power_flow.parameters = {
                'numBuses', 'busData', 'lineData', 'method', 'tolerance', 'maxIter'
            };
            domains.power_flow.subCategories = {
                'Load Flow', 'Optimal Power Flow', 'DC Power Flow', 'Continuation Power Flow'
            };
            domains.power_flow.methods = {
                'Newton-Raphson', 'Gauss-Seidel', 'Fast Decoupled', 'DC Power Flow'
            };

            % Fault Analysis
            domains.fault_analysis.name = 'Fault Analysis';
            domains.fault_analysis.description = 'Short circuit studies and fault calculations';
            domains.fault_analysis.module = 'FaultAnalysisModule';
            domains.fault_analysis.implementation = 'code';
            domains.fault_analysis.parameters = {
                'faultType', 'faultLocation', 'faultImpedance', 'XR_ratio', 'voltageLevel'
            };
            domains.fault_analysis.subCategories = {
                'Three-Phase Fault', 'Single Line-to-Ground', 'Line-to-Line', 'Double Line-to-Ground'
            };
            domains.fault_analysis.methods = {
                'Impedance Matrix', 'Symmetrical Components', 'Sequence Networks'
            };

            % Protection Systems
            domains.protection.name = 'Protection Systems';
            domains.protection.description = 'Relay coordination and protection design';
            domains.protection.module = 'ProtectionModule';
            domains.protection.implementation = 'simulink';
            domains.protection.parameters = {
                'relayType', 'CT_ratio', 'VT_ratio', 'pickupCurrent', 'TDS', 'CTI'
            };
            domains.protection.subCategories = {
                'Overcurrent Protection', 'Distance Protection', 'Differential Protection', 'Relay Coordination'
            };
            domains.protection.methods = {
                'Time-Current Coordination', 'Zone Coordination', 'TCC Analysis'
            };

            % Stability Studies
            domains.stability.name = 'Stability Studies';
            domains.stability.description = 'System stability analysis and assessment';
            domains.stability.module = 'StabilityModule';
            domains.stability.implementation = 'simulink';
            domains.stability.parameters = {
                'studyType', 'simulationTime', 'stepSize', 'generatorData', 'disturbanceType'
            };
            domains.stability.subCategories = {
                'Transient Stability', 'Small Signal Stability', 'Voltage Stability', 'Frequency Stability'
            };
            domains.stability.methods = {
                'Time-Domain Simulation', 'Eigenvalue Analysis', 'P-V/Q-V Curves'
            };

            % Renewable Integration
            domains.renewable_integration.name = 'Renewable Integration';
            domains.renewable_integration.description = 'Modeling of renewable energy sources';
            domains.renewable_integration.module = 'RenewableModule';
            domains.renewable_integration.implementation = 'simulink';
            domains.renewable_integration.parameters = {
                'resourceType', 'capacity', 'penetrationLevel', 'controlMode', 'storageCapacity'
            };
            domains.renewable_integration.subCategories = {
                'Solar PV Systems', 'Wind Farms', 'Energy Storage', 'Microgrids', 'Grid-Tied Inverters'
            };
            domains.renewable_integration.methods = {
                'Phasor Modeling', 'Electromagnetic Transient', 'Average Value Modeling'
            };
        end

        function category = classifyByTask(obj, taskDescription)
            % Classify a task into appropriate power system category
            taskLower = lower(taskDescription);

            % Classification keywords
            categories = {
                'power_flow', 'fault_analysis', 'protection', 'stability', 'renewable_integration'
            };

            keywords = containers.Map();
            keywords('power_flow') = {
                'load flow', 'power flow', 'voltage profile', 'power transfer', ...
                'bus voltage', 'line loading', 'reactive power', 'PV bus', 'PQ bus'
            };
            keywords('fault_analysis') = {
                'fault', 'short circuit', 'fault current', 'symmetrical components', ...
                'sequence impedance', 'fault level', 'breaker rating'
            };
            keywords('protection') = {
                'relay', 'protection', 'coordination', 'overcurrent', 'distance relay', ...
                'differential', 'TCC', 'trip', 'CT', 'VT', 'pickup'
            };
            keywords('stability') = {
                'stability', 'oscillation', 'damping', 'rotor angle', 'critical clearing', ...
                'eigenvalue', 'small signal', 'transient stability', 'voltage collapse'
            };
            keywords('renewable_integration') = {
                'renewable', 'solar', 'wind', 'PV', 'inverter', 'distributed generation', ...
                'microgrid', 'energy storage', 'battery', 'DFIG', 'grid-tied'
            };

            % Score each category
            scores = zeros(1, length(categories));
            for i = 1:length(categories)
                catKeywords = keywords(categories{i});
                for j = 1:length(catKeywords)
                    if contains(taskLower, catKeywords{j})
                        scores(i) = scores(i) + 1;
                    end
                end
            end

            % Find best match
            [maxScore, maxIdx] = max(scores);
            if maxScore > 0
                category = categories{maxIdx};
            else
                category = 'unknown';
            end
        end

        function suggestions = getParameterSuggestions(obj, domain)
            % Get suggested parameters for a domain
            suggestions = struct();

            switch domain
                case 'power_flow'
                    suggestions.method = 'Newton-Raphson';
                    suggestions.tolerance = 1e-6;
                    suggestions.maxIter = 50;
                    suggestions.baseMVA = 100;

                case 'fault_analysis'
                    suggestions.faultType = 'three_phase';
                    suggestions.faultImpedance = 0.01;
                    suggestions.XR_ratio = 10;
                    suggestions.baseMVA = 100;

                case 'protection'
                    suggestions.relayType = 'overcurrent';
                    suggestions.CTI = 0.3;
                    suggestions.TDS = 1.0;
                    suggestions.curveType = 'IEEE moderately inverse';

                case 'stability'
                    suggestions.studyType = 'transient';
                    suggestions.simulationTime = 10;
                    suggestions.stepSize = 0.001;
                    suggestions.solver = 'ode23t';

                case 'renewable_integration'
                    suggestions.resourceType = 'solar_pv';
                    suggestions.controlMode = 'PQ';
                    suggestions.penetrationLevel = 20;
                    suggestions.includeStorage = false;

                otherwise
                    suggestions = [];
            end
        end

        function validValues = getValidValues(obj, paramName)
            % Get valid values/options for a parameter
            validOptions = obj.getValidParameterOptions();

            if isfield(validOptions, paramName)
                validValues = validOptions.(paramName);
            else
                validValues = [];
            end
        end

        function options = getValidParameterOptions(obj)
            % Get all valid parameter options
            options.method = {'Newton-Raphson', 'Gauss-Seidel', 'Fast Decoupled', 'DC Power Flow'};
            options.faultType = {'three_phase', 'single_line_ground', 'line_to_line', 'double_line_ground'};
            options.relayType = {'overcurrent', 'distance', 'differential', 'directional'};
            options.studyType = {'transient', 'small_signal', 'voltage', 'frequency'};
            options.resourceType = {'solar_pv', 'wind', 'storage', 'hybrid'};
            options.controlMode = {'PQ', 'PV', 'droop', 'grid_forming'};
            options.solver = {'ode45', 'ode23t', 'ode15s', 'ode23'};
            options.curveType = {
                'IEEE moderately inverse', 'IEEE very inverse', 'IEEE extremely inverse', ...
                'IEC normal inverse', 'IEC very inverse', 'IEC extremely inverse'
            };
        end

        function relatedDomains = getRelatedDomains(obj, domain)
            % Get related domains for cross-domain analysis
            relations = containers.Map();
            relations('power_flow') = {'stability', 'renewable_integration'};
            relations('fault_analysis') = {'protection', 'stability'};
            relations('protection') = {'fault_analysis', 'stability'};
            relations('stability') = {'power_flow', 'renewable_integration', 'protection'};
            relations('renewable_integration') = {'power_flow', 'stability', 'control_systems'};

            if isKey(relations, domain)
                relatedDomains = relations(domain);
            else
                relatedDomains = {};
            end
        end
    end
end
