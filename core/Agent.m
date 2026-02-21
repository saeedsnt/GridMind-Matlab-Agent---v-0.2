classdef Agent < handle
    % Agent - Main controller for the MATLAB Power Systems AI Agent
    %
    % This class orchestrates the interactive session, handling user input,
    % classification, and coordination of domain-specific modules.

    properties
        currentDomain           % Current power system domain
        implementationType      % 'code' or 'simulink'
        parameters              % Collected simulation parameters
        recommendationEngine    % Best practices engine
        interactivePrompt       % User input handler
        activeModule            % Current domain module
    end

    methods
        function obj = Agent()
            % Constructor - initialize components
            obj.recommendationEngine = RecommendationEngine();
            obj.interactivePrompt = InteractivePrompt();
            obj.parameters = struct();
            obj.currentDomain = [];
            obj.implementationType = [];
            obj.activeModule = [];
        end

        function startInteractiveSession(obj)
            % Start the main interactive loop
            while true
                try
                    % Get user input
                    userInput = obj.interactivePrompt.getInput();

                    % Process input
                    if obj.processInput(userInput)
                        break;  % Exit on quit command
                    end

                catch ME
                    fprintf('\nError: %s\n', ME.message);
                    fprintf('Please try again or type /help for assistance.\n\n');
                end
            end
        end

        function shouldExit = processInput(obj, userInput)
            % Process user input and return true if should exit
            shouldExit = false;

            % Trim and convert to lowercase for comparison
            cmd = strtrim(lower(userInput));

            % Handle slash commands
            if startsWith(cmd, '/')
                switch cmd
                    case '/quit'
                        obj.displayGoodbye();
                        shouldExit = true;
                        return;

                    case '/help'
                        obj.displayHelp();
                        return;

                    case '/power-flow'
                        obj.initializeDomain('power_flow');

                    case '/power-fault'
                        obj.initializeDomain('fault_analysis');

                    case '/power-stability'
                        obj.initializeDomain('stability');

                    case '/power-protection'
                        obj.initializeDomain('protection');

                    case '/power-renewable'
                        obj.initializeDomain('renewable_integration');

                    otherwise
                        fprintf('Unknown command: %s\n', cmd);
                        fprintf('Type /help for available commands.\n\n');
                end
            else
                % Handle contextual input based on current state
                obj.handleContextualInput(userInput);
            end
        end

        function initializeDomain(obj, domainName)
            % Initialize a specific power system domain
            fprintf('\nInitializing %s domain...\n', strrep(domainName, '_', ' '));

            % Create appropriate classifier
            classifier = PowerSystemClassifier();
            domainInfo = classifier.getDomainInfo(domainName);

            if isempty(domainInfo)
                fprintf('Domain not found: %s\n', domainName);
                return;
            end

            obj.currentDomain = domainName;

            % Determine implementation type
            implClassifier = ImplementationClassifier();
            obj.implementationType = implClassifier.suggestImplementation(domainName);

            % Ask user for implementation preference
            obj.promptImplementationChoice();

            % Create domain module
            obj.createDomainModule(domainName);

            % Start parameter collection
            obj.collectParameters();
        end

        function promptImplementationChoice(obj)
            % Let user choose or confirm implementation type
            fprintf('\n--------------------------------------------------------------------------------\n');
            fprintf('Implementation Type Selection\n');
            fprintf('--------------------------------------------------------------------------------\n');
            fprintf('Based on your selection, %s is recommended.\n', obj.implementationType);
            fprintf('\n');
            fprintf('  1. Code-based (MATLAB scripts for calculations and analysis)\n');
            fprintf('  2. Simulink-based (Block diagram simulations and dynamic modeling)\n');
            fprintf('\n');

            choice = obj.interactivePrompt.getNumericInput('Choose implementation (1-2): ', [1, 2]);

            if choice == 1
                obj.implementationType = 'code';
            else
                obj.implementationType = 'simulink';
            end

            fprintf('Selected: %s-based implementation\n\n', obj.implementationType);
        end

        function createDomainModule(obj, domainName)
            % Create the appropriate domain module
            switch domainName
                case 'power_flow'
                    obj.activeModule = PowerFlowModule();
                case 'fault_analysis'
                    obj.activeModule = FaultAnalysisModule();
                case 'protection'
                    obj.activeModule = ProtectionModule();
                case 'stability'
                    obj.activeModule = StabilityModule();
                case 'renewable_integration'
                    obj.activeModule = RenewableModule();
                otherwise
                    error('Unknown domain: %s', domainName);
            end
        end

        function collectParameters(obj)
            % Collect parameters from user with guidance
            if isempty(obj.activeModule)
                error('No active module for parameter collection');
            end

            fprintf('\n--------------------------------------------------------------------------------\n');
            fprintf('Parameter Collection for %s\n', strrep(obj.currentDomain, '_', ' '));
            fprintf('--------------------------------------------------------------------------------\n\n');

            % Get module-specific parameters
            obj.parameters = obj.activeModule.collectParameters(obj.interactivePrompt);

            % Validate parameters
            validator = Validator();
            validationResults = validator.validateParameters(obj.parameters, obj.currentDomain);

            if ~validationResults.isValid
                fprintf('\nValidation warnings:\n');
                for i = 1:length(validationResults.warnings)
                    fprintf('  - %s\n', validationResults.warnings{i});
                end

                retry = obj.interactivePrompt.getYesNo('Would you like to modify parameters? (y/n): ');
                if retry
                    obj.collectParameters();
                    return;
                end
            end

            % Generate code
            obj.generateOutput();
        end

        function generateOutput(obj)
            % Generate MATLAB/Simulink code based on collected parameters
            fprintf('\n--------------------------------------------------------------------------------\n');
            fprintf('Generating Code\n');
            fprintf('--------------------------------------------------------------------------------\n\n');

            if isempty(obj.activeModule)
                error('No active module for code generation');
            end

            % Generate code using active module
            if strcmp(obj.implementationType, 'code')
                generatedCode = obj.activeModule.generateCode(obj.parameters);

                % Clean code of any branding
                cleaner = CodeCleaner();
                generatedCode = cleaner.clean(generatedCode);

                % Display generated code
                fprintf('Generated MATLAB Code:\n');
                fprintf('================================================================================\n');
                fprintf('%s\n', generatedCode);
                fprintf('================================================================================\n\n');

                % Offer to save
                saveFile = obj.interactivePrompt.getYesNo('Save to file? (y/n): ');
                if saveFile
                    defaultName = [obj.currentDomain '_script.m'];
                    fileName = obj.interactivePrompt.getInput(['Enter filename (default: ' defaultName '): ']);
                    if isempty(fileName)
                        fileName = defaultName;
                    end
                    if ~endsWith(fileName, '.m')
                        fileName = [fileName '.m'];
                    end

                    fid = fopen(fileName, 'w');
                    fprintf(fid, '%s', generatedCode);
                    fclose(fid);
                    fprintf('Code saved to: %s\n\n', fileName);
                end

            else
                % Simulink model generation
                fprintf('Simulink model generation:\n');
                obj.activeModule.generateSimulinkModel(obj.parameters);
            end

            % Offer recommendations
            obj.showRecommendations();

            % Ask if user wants to continue
            continueSession = obj.interactivePrompt.getYesNo('Start another analysis? (y/n): ');
            if ~continueSession
                obj.displayGoodbye();
                exit;
            end
        end

        function showRecommendations(obj)
            % Show best practice recommendations
            fprintf('\n--------------------------------------------------------------------------------\n');
            fprintf('Best Practice Recommendations\n');
            fprintf('--------------------------------------------------------------------------------\n\n');

            recommendations = obj.recommendationEngine.getRecommendations(obj.currentDomain, obj.parameters);

            for i = 1:length(recommendations)
                fprintf('  [Tip %d] %s\n', i, recommendations{i});
            end
            fprintf('\n');
        end

        function handleContextualInput(obj, userInput)
            % Handle input based on current context
            if isempty(obj.currentDomain)
                fprintf('Please select a domain using one of the slash commands.\n');
                fprintf('Type /help for available commands.\n\n');
            else
                % Context-specific handling delegated to active module
                if ~isempty(obj.activeModule)
                    obj.activeModule.handleInput(userInput);
                end
            end
        end

        function displayHelp(obj)
            % Display help information
            fprintf('\n');
            fprintf('================================================================================\n');
            fprintf('                           Help - Available Commands\n');
            fprintf('================================================================================\n');
            fprintf('\n');
            fprintf('Domain Commands:\n');
            fprintf('  /power-flow       - Configure and run load flow analysis\n');
            fprintf('  /power-fault      - Set up short circuit and fault studies\n');
            fprintf('  /power-stability  - Analyze system stability margins\n');
            fprintf('  /power-protection - Design and coordinate protection relays\n');
            fprintf('  /power-renewable  - Model renewable energy integration\n');
            fprintf('\n');
            fprintf('General Commands:\n');
            fprintf('  /help             - Display this help message\n');
            fprintf('  /quit             - Exit the agent\n');
            fprintf('\n');
            fprintf('Tips:\n');
            fprintf('  - The agent will guide you through parameter collection\n');
            fprintf('  - Default values follow IEEE and IEC standards\n');
            fprintf('  - You can modify parameters before final code generation\n');
            fprintf('\n');
            fprintf('================================================================================\n');
            fprintf('\n');
        end

        function displayGoodbye(obj)
            % Display goodbye message
            fprintf('\n');
            fprintf('--------------------------------------------------------------------------------\n');
            fprintf('  Thank you for using the MATLAB AI Agent for Power Systems.\n');
            fprintf('  Goodbye!\n');
            fprintf('--------------------------------------------------------------------------------\n');
            fprintf('\n');
        end
    end
end
