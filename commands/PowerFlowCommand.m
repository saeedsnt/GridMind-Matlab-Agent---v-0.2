classdef PowerFlowCommand < handle
    % PowerFlowCommand - Handler for power flow analysis commands
    %
    % Manages power flow analysis workflows including parameter setup,
    % validation, and execution.

    properties
        module              % PowerFlowModule instance
        parameters          % Collected parameters
        results             % Analysis results
    end

    methods
        function obj = PowerFlowCommand()
            % Constructor
            obj.module = PowerFlowModule();
            obj.parameters = struct();
            obj.results = struct();
        end

        function execute(obj, prompt)
            % Execute power flow analysis
            fprintf('\n--- Power Flow Analysis Command ---\n');
            obj.parameters = obj.module.collectParameters(prompt);
            obj.displayResults();
        end

        function displayResults(obj)
            % Display analysis results
            fprintf('\n=== Analysis Configuration Complete ===\n');
            fprintf('Method: %s\n', obj.parameters.method);
            fprintf('Buses: %d\n', obj.parameters.numBuses);
            fprintf('Base MVA: %.2f\n', obj.parameters.baseMVA);
            fprintf('Tolerance: %g\n', obj.parameters.tolerance);
            fprintf('\n');
        end
    end
end
