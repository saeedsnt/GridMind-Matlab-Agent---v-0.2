%% MATLAB AI Agent for Power Systems
% Main entry point - Interactive power system simulation setup
%
% This agent provides an intelligent, step-by-step simulation setup
% experience with best practice recommendations for power system engineers.

function main()
    % Initialize the agent
    agent = Agent();

    % Display welcome message
    displayWelcome();

    % Start interactive session
    agent.startInteractiveSession();
end

function displayWelcome()
    % Clear command window for clean display
    clc;

    % Display banner
    fprintf('\n');
    fprintf('================================================================================\n');
    fprintf('           MATLAB AI Agent for Power Systems\n');
    fprintf('================================================================================\n');
    fprintf('\n');
    fprintf('  Welcome to the intelligent power system simulation assistant.\n');
    fprintf('  This tool guides you through simulation setup with best practice\n');
    fprintf('  recommendations and industry-standard parameters.\n');
    fprintf('\n');
    fprintf('--------------------------------------------------------------------------------\n');
    fprintf('  Available Commands:\n');
    fprintf('--------------------------------------------------------------------------------\n');
    fprintf('  /power-flow       - Load flow and power flow analysis\n');
    fprintf('  /power-fault      - Short circuit and fault analysis\n');
    fprintf('  /power-stability  - System stability studies\n');
    fprintf('  /power-protection - Relay coordination and protection\n');
    fprintf('  /power-renewable  - Renewable energy integration\n');
    fprintf('  /help             - Show this help message\n');
    fprintf('  /quit             - Exit the agent\n');
    fprintf('--------------------------------------------------------------------------------\n');
    fprintf('\n');
end
