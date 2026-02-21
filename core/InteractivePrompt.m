classdef InteractivePrompt < handle
    % InteractivePrompt - Handles user input for the power systems agent
    %
    % Provides methods for getting various types of user input with
    % validation and user-friendly prompts.

    methods
        function obj = InteractivePrompt()
            % Constructor
        end

        function userInput = getInput(obj, promptText)
            % Get string input from user
            if nargin < 2
                promptText = '> ';
            end
            userInput = input(promptText, 's');
        end

        function numericValue = getNumericInput(obj, promptText, validRange)
            % Get numeric input with optional range validation
            if nargin < 2
                promptText = 'Enter value: ';
            end

            while true
                try
                    str = input(promptText, 's');
                    numericValue = str2double(str);

                    if isnan(numericValue)
                        fprintf('Please enter a valid number.\n');
                        continue;
                    end

                    % Check range if provided
                    if nargin >= 3 && ~isempty(validRange)
                        if length(validRange) == 2
                            if numericValue < validRange(1) || numericValue > validRange(2)
                                fprintf('Value must be between %.2f and %.2f.\n', ...
                                    validRange(1), validRange(2));
                                continue;
                            end
                        elseif length(validRange) == 1
                            if numericValue ~= validRange(1)
                                fprintf('Value must be %.2f.\n', validRange(1));
                                continue;
                            end
                        end
                    end

                    return;

                catch ME
                    fprintf('Invalid input. Please try again.\n');
                end
            end
        end

        function numericValue = getNumericInputWithDefault(obj, promptText, defaultValue, validRange)
            % Get numeric input with default value and optional range validation
            if nargin < 2
                promptText = 'Enter value: ';
            end

            defaultStr = sprintf('%.4g', defaultValue);
            fullPrompt = sprintf('%s [default: %s]: ', promptText, defaultStr);

            while true
                try
                    str = input(fullPrompt, 's');

                    % Use default if empty
                    if isempty(strtrim(str))
                        numericValue = defaultValue;
                        return;
                    end

                    numericValue = str2double(str);

                    if isnan(numericValue)
                        fprintf('Please enter a valid number.\n');
                        continue;
                    end

                    % Check range if provided
                    if nargin >= 4 && ~isempty(validRange)
                        if length(validRange) == 2
                            if numericValue < validRange(1) || numericValue > validRange(2)
                                fprintf('Value must be between %.2f and %.2f.\n', ...
                                    validRange(1), validRange(2));
                                continue;
                            end
                        end
                    end

                    return;

                catch ME
                    fprintf('Invalid input. Please try again.\n');
                end
            end
        end

        function yesNo = getYesNo(obj, promptText)
            % Get yes/no input from user
            if nargin < 2
                promptText = 'Continue? (y/n): ';
            end

            while true
                str = input(promptText, 's');
                str = lower(strtrim(str));

                if isempty(str)
                    yesNo = true;  % Default to yes
                    return;
                end

                if startsWith(str, 'y')
                    yesNo = true;
                    return;
                elseif startsWith(str, 'n')
                    yesNo = false;
                    return;
                else
                    fprintf('Please enter y or n.\n');
                end
            end
        end

        function choice = getMenuChoice(obj, options, promptText)
            % Get menu choice from list of options
            if nargin < 3
                promptText = 'Select option: ';
            end

            % Display options
            fprintf('\n');
            for i = 1:length(options)
                fprintf('  [%d] %s\n', i, options{i});
            end
            fprintf('\n');

            while true
                choice = obj.getNumericInput(promptText, [1, length(options)]);
                return;
            end
        end

        function choice = getMenuChoiceWithDefault(obj, options, defaultIdx, promptText)
            % Get menu choice with default option
            if nargin < 4
                promptText = 'Select option: ';
            end

            % Display options with default indicator
            fprintf('\n');
            for i = 1:length(options)
                if i == defaultIdx
                    fprintf('  [%d] %s (default)\n', i, options{i});
                else
                    fprintf('  [%d] %s\n', i, options{i});
                end
            end
            fprintf('\n');

            defaultStr = sprintf('[default: %d] ', defaultIdx);
            fullPrompt = [defaultStr promptText];

            while true
                str = input(fullPrompt, 's');

                if isempty(strtrim(str))
                    choice = defaultIdx;
                    return;
                end

                choice = str2double(str);
                if isnan(choice) || choice < 1 || choice > length(options)
                    fprintf('Please enter a number between 1 and %d.\n', length(options));
                    continue;
                end
                return;
            end
        end

        function values = getVectorInput(obj, promptText, size, defaultValue)
            % Get vector/matrix input from user
            if nargin < 2
                promptText = 'Enter values: ';
            end

            if nargin >= 4 && ~isempty(defaultValue)
                defStr = mat2str(defaultValue);
                fullPrompt = sprintf('%s [default: %s]: ', promptText, defStr);
            else
                fullPrompt = promptText;
            end

            while true
                str = input(fullPrompt, 's');

                if isempty(strtrim(str)) && nargin >= 4 && ~isempty(defaultValue)
                    values = defaultValue;
                    return;
                end

                try
                    values = str2num(str);
                    if isempty(values)
                        fprintf('Invalid input. Please try again.\n');
                        continue;
                    end
                    return;
                catch
                    fprintf('Invalid input. Please try again.\n');
                end
            end
        end

        function value = getSelectOne(obj, options, descriptions, promptText)
            % Get single selection from options with descriptions
            if nargin < 4
                promptText = 'Select option: ';
            end

            fprintf('\n');
            for i = 1:length(options)
                fprintf('  [%d] %s\n', i, options{i});
                if nargin >= 3 && ~isempty(descriptions)
                    fprintf('      %s\n', descriptions{i});
                end
            end
            fprintf('\n');

            while true
                choice = input([promptText ' [1-' num2str(length(options)) ']: '], 's');
                choice = str2double(choice);

                if isnan(choice) || choice < 1 || choice > length(options)
                    fprintf('Invalid selection. Please try again.\n');
                    continue;
                end

                value = options{choice};
                return;
            end
        end
    end
end
