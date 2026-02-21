classdef RenewableIntegrationModule < handle
    % RenewableIntegrationModule - Renewable energy integration domain module
    %
    % Provides interactive parameter collection and code generation
    % for modeling and analysis of renewable energy integration including
    % solar PV, wind, storage, and microgrid systems.

    properties
        resourceType         % Type of renewable resource
        capacity             % System capacity (MW)
        penetrationLevel     % Penetration level (%)
        controlMode          % Control strategy
        storageCapacity      % Energy storage capacity (MWh)
        inverterType         % Type of inverter
        gridFormingEnabled   % Grid-forming capability
        siteLongitude        % Site longitude for solar modeling
        siteLatitude         % Site latitude for solar modeling
        weatherData          % Wind/solar weather data
        systemVoltage        % System voltage level (kV)
        baseMVA              % Base MVA for per-unit system
    end

    methods
        function obj = RenewableIntegrationModule()
            % Constructor - set defaults
            obj.resourceType = 'solar_pv';
            obj.capacity = 10;
            obj.penetrationLevel = 20;
            obj.controlMode = 'PQ';
            obj.storageCapacity = 5;
            obj.inverterType = 'string_inverter';
            obj.gridFormingEnabled = false;
            obj.siteLongitude = -80;
            obj.siteLatitude = 40;
            obj.weatherData = [];
            obj.systemVoltage = 13.8;
            obj.baseMVA = 100;
        end

        function params = collectParameters(obj, prompt)
            % Collect renewable integration parameters from user

            fprintf('\n=== Renewable Integration Configuration ===\n\n');

            % Select resource type
            resourceTypes = {'solar_pv', 'wind', 'storage', 'hybrid_pv_wind', 'microgrid'};
            fprintf('Available renewable resource types:\n');
            resourceIdx = prompt.getMenuChoice(resourceTypes, 'Select resource type:');
            obj.resourceType = resourceTypes{resourceIdx};
            fprintf('Selected resource type: %s\n\n', obj.resourceType);

            % System capacity
            obj.capacity = prompt.getNumericInputWithDefault( ...
                'Enter system capacity (MW): ', 10, [0.1, 1000]);

            % Penetration level
            obj.penetrationLevel = prompt.getNumericInputWithDefault( ...
                'Enter penetration level (%): ', 20, [1, 100]);

            % Control mode
            controlModes = {'PQ', 'PV', 'droop', 'grid_forming'};
            controlIdx = prompt.getMenuChoice(controlModes, 'Select control mode:');
            obj.controlMode = controlModes{controlIdx};
            fprintf('Selected control mode: %s\n\n', obj.controlMode);

            % Inverter type based on resource
            switch obj.resourceType
                case 'solar_pv'
                    fprintf('\n--- Solar PV Configuration ---\n');
                    obj = obj.collectSolarParameters(prompt);

                case 'wind'
                    fprintf('\n--- Wind Farm Configuration ---\n');
                    obj = obj.collectWindParameters(prompt);

                case 'storage'
                    fprintf('\n--- Battery Energy Storage Configuration ---\n');
                    obj = obj.collectStorageParameters(prompt);

                case 'hybrid_pv_wind'
                    fprintf('\n--- Hybrid PV+Wind Configuration ---\n');
                    obj = obj.collectHybridParameters(prompt);

                case 'microgrid'
                    fprintf('\n--- Microgrid Configuration ---\n');
                    obj = obj.collectMicrogridParameters(prompt);
            end

            % System and grid parameters
            fprintf('\n--- System Parameters ---\n');
            obj.systemVoltage = prompt.getNumericInputWithDefault( ...
                'Enter system voltage (kV): ', 13.8, [0.208, 765]);

            obj.baseMVA = prompt.getNumericInputWithDefault( ...
                'Enter base MVA: ', 100, [1, 10000]);

            % Grid-forming capability
            enableGridForming = prompt.getYesNo('Enable grid-forming capability? (y/n): ');
            obj.gridFormingEnabled = enableGridForming;

            % Return parameters struct
            params.resourceType = obj.resourceType;
            params.capacity = obj.capacity;
            params.penetrationLevel = obj.penetrationLevel;
            params.controlMode = obj.controlMode;
            params.storageCapacity = obj.storageCapacity;
            params.inverterType = obj.inverterType;
            params.gridFormingEnabled = obj.gridFormingEnabled;
            params.systemVoltage = obj.systemVoltage;
            params.baseMVA = obj.baseMVA;
        end

        function obj = collectSolarParameters(obj, prompt)
            % Collect solar PV specific parameters
            obj.inverterType = 'string_inverter';

            % Array configuration
            numModules = prompt.getNumericInputWithDefault( ...
                'Enter number of PV modules: ', 1000, [1, 100000]);

            moduleRating = prompt.getNumericInputWithDefault( ...
                'Enter module rating (W): ', 350, [100, 1000]);

            % Temperature coefficient
            tempCoeff = prompt.getNumericInputWithDefault( ...
                'Enter temperature coefficient (%/°C): ', -0.4, [-0.6, -0.2]);

            % Reference conditions
            STC = prompt.getNumericInputWithDefault( ...
                'Reference irradiance (W/m²): ', 1000, [100, 1200]);

            fprintf('Solar system configured: %.2f kW installed capacity\n\n', ...
                numModules * moduleRating / 1000);

            % Store in struct for code generation
            obj.weatherData = struct(...
                'numModules', numModules, ...
                'moduleRating', moduleRating, ...
                'tempCoeff', tempCoeff, ...
                'STC', STC);
        end

        function obj = collectWindParameters(obj, prompt)
            % Collect wind farm specific parameters
            generatorTypes = {'DFIG', 'PMSG', 'SG'};
            fprintf('Wind generator types: 1=DFIG, 2=PMSG, 3=SG\n');
            genTypeIdx = prompt.getMenuChoice(generatorTypes, 'Select wind generator type:');
            obj.inverterType = generatorTypes{genTypeIdx};

            % Turbine parameters
            numTurbines = prompt.getNumericInputWithDefault( ...
                'Enter number of turbines: ', 10, [1, 500]);

            ratedPower = prompt.getNumericInputWithDefault( ...
                'Enter turbine rated power (MW): ', 2.5, [0.5, 15]);

            % Wind characteristics
            hubHeight = prompt.getNumericInputWithDefault( ...
                'Enter hub height (m): ', 90, [20, 150]);

            rotorDiameter = prompt.getNumericInputWithDefault( ...
                'Enter rotor diameter (m): ', 100, [30, 200]);

            % Wind speed parameters
            minWindSpeed = prompt.getNumericInputWithDefault( ...
                'Enter minimum wind speed (m/s): ', 3, [1, 5]);

            maxWindSpeed = prompt.getNumericInputWithDefault( ...
                'Enter maximum wind speed (m/s): ', 25, [15, 35]);

            fprintf('Wind farm configured: %.2f MW installed capacity\n\n', ...
                numTurbines * ratedPower);

            obj.weatherData = struct(...
                'numTurbines', numTurbines, ...
                'ratedPower', ratedPower, ...
                'hubHeight', hubHeight, ...
                'rotorDiameter', rotorDiameter, ...
                'minWindSpeed', minWindSpeed, ...
                'maxWindSpeed', maxWindSpeed);
        end

        function obj = collectStorageParameters(obj, prompt)
            % Collect battery energy storage specific parameters
            obj.inverterType = 'bidirectional_converter';

            % Battery specifications
            obj.storageCapacity = prompt.getNumericInputWithDefault( ...
                'Enter energy storage capacity (MWh): ', 5, [0.1, 1000]);

            nominalPower = prompt.getNumericInputWithDefault( ...
                'Enter nominal power (MW): ', obj.capacity, [0.1, 1000]);

            % Battery chemistry
            chemistries = {'lithium_ion', 'lead_acid', 'flow_battery', 'supercapacitor'};
            chemistrySel = prompt.getMenuChoice(chemistries, 'Select battery chemistry:');
            chemistry = chemistries{chemistrySel};

            % Efficiency parameters
            roundTripEfficiency = prompt.getNumericInputWithDefault( ...
                'Enter round-trip efficiency (%): ', 85, [60, 95]);

            chargeTime = prompt.getNumericInputWithDefault( ...
                'Enter charge time (hours): ', 2, [0.5, 24]);

            dischargeTime = prompt.getNumericInputWithDefault( ...
                'Enter discharge time (hours): ', 2, [0.5, 24]);

            fprintf('Battery storage configured: %.2f MWh capacity\n\n', obj.storageCapacity);

            obj.weatherData = struct(...
                'nominalPower', nominalPower, ...
                'chemistry', chemistry, ...
                'roundTripEfficiency', roundTripEfficiency, ...
                'chargeTime', chargeTime, ...
                'dischargeTime', dischargeTime, ...
                'stateOfCharge', 50);  % Initial SOC
        end

        function obj = collectHybridParameters(obj, prompt)
            % Collect hybrid PV+Wind system parameters
            fprintf('Configuring hybrid PV+Wind system...\n\n');

            % Solar portion
            fprintf('--- Solar PV Configuration ---\n');
            solarCapacity = prompt.getNumericInputWithDefault( ...
                'Enter solar PV capacity (MW): ', obj.capacity/2, [0.1, 500]);

            % Wind portion
            fprintf('\n--- Wind Farm Configuration ---\n');
            windCapacity = prompt.getNumericInputWithDefault( ...
                'Enter wind farm capacity (MW): ', obj.capacity/2, [0.1, 500]);

            % Storage
            fprintf('\n--- Energy Storage Configuration ---\n');
            obj.storageCapacity = prompt.getNumericInputWithDefault( ...
                'Enter energy storage capacity (MWh): ', 5, [0.1, 500]);

            obj.inverterType = 'hybrid_inverter';

            fprintf('\nHybrid system configured: %.2f MW PV, %.2f MW Wind, %.2f MWh Storage\n\n', ...
                solarCapacity, windCapacity, obj.storageCapacity);

            obj.weatherData = struct(...
                'solarCapacity', solarCapacity, ...
                'windCapacity', windCapacity, ...
                'storageCapacity', obj.storageCapacity);
        end

        function obj = collectMicrogridParameters(obj, prompt)
            % Collect microgrid system parameters
            fprintf('Configuring microgrid system...\n\n');

            % Number of DG units
            numDGs = prompt.getNumericInputWithDefault( ...
                'Enter number of DG units: ', 3, [1, 20]);

            % Control architecture
            controlArchs = {'centralized', 'decentralized', 'hierarchical'};
            controlIdx = prompt.getMenuChoice(controlArchs, 'Select control architecture:');
            controlArch = controlArchs{controlIdx};

            % Island capability
            canIsland = prompt.getYesNo('Enable island mode capability? (y/n): ');

            % Loads
            totalLoad = prompt.getNumericInputWithDefault( ...
                'Enter total microgrid load (MW): ', 5, [0.1, 100]);

            % Communication latency
            commLatency = prompt.getNumericInputWithDefault( ...
                'Enter communication latency (ms): ', 10, [1, 1000]);

            fprintf('Microgrid configured: %d DG units, %.2f MW load\n\n', numDGs, totalLoad);

            obj.weatherData = struct(...
                'numDGs', numDGs, ...
                'controlArchitecture', controlArch, ...
                'canIsland', canIsland, ...
                'totalLoad', totalLoad, ...
                'commLatency', commLatency);
        end

        function code = generateCode(obj, params)
            % Generate MATLAB code for renewable integration analysis
            code = obj.generateRenewableScript();
        end

        function script = generateRenewableScript(obj)
            % Generate complete renewable integration script

            script = [
                '%% Renewable Energy Integration Analysis' newline
                '% Auto-generated MATLAB code for renewable energy system analysis' newline
                newline
                'clear; clc; close all;' newline
                newline
                '%% System Parameters' newline
                sprintf('capacity = %.2f;  %% MW' newline, obj.capacity)
                sprintf('penetrationLevel = %.2f;  %% %%' newline, obj.penetrationLevel)
                sprintf('controlMode = ''%s'';' newline, obj.controlMode)
                sprintf('resourceType = ''%s'';' newline, obj.resourceType)
                sprintf('systemVoltage = %.2f;  %% kV' newline, obj.systemVoltage)
                sprintf('baseMVA = %.2f;' newline, obj.baseMVA)
                newline
            ];

            % Add resource-specific code
            switch obj.resourceType
                case 'solar_pv'
                    script = [script, obj.getSolarCode()];
                case 'wind'
                    script = [script, obj.getWindCode()];
                case 'storage'
                    script = [script, obj.getStorageCode()];
                case 'hybrid_pv_wind'
                    script = [script, obj.getHybridCode()];
                case 'microgrid'
                    script = [script, obj.getMicrogridCode()];
            end

            % Add common analysis code
            script = [script, obj.getCommonAnalysisCode()];
        end

        function code = getSolarCode(obj)
            % Generate solar PV specific code
            code = [
                '%% Solar PV System Model' newline
                newline
                'numModules = %d;' newline, obj.weatherData.numModules
                'moduleRating = %d;  %% W' newline, obj.weatherData.moduleRating
                'tempCoeff = %.4f;  %% per degree C' newline, obj.weatherData.tempCoeff
                'STC_irradiance = %.0f;  %% W/m²' newline, obj.weatherData.STC
                'T_ref = 25;  %% Reference temperature (°C)' newline
                newline
                '%% Irradiance Profile (Example)' newline
                '% Typical solar irradiance on clear day' newline
                'hours = 0:24;' newline
                'irradiance = 1000 * max(0, sin(pi * hours / 24)).^1.5;' newline
                'temperature = 25 + 10 * sin(pi * hours / 24);  %% Temperature variation' newline
                newline
                '%% PV Array Output Calculation' newline
                'V_oc = 37.2;  %% Open-circuit voltage (V) per module' newline
                'I_sc = 8.5;   %% Short-circuit current (A)' newline
                'V_mp = 30.5;  %% Max power point voltage' newline
                'I_mp = 8.0;   %% Max power point current' newline
                'P_mp = V_mp * I_mp;  %% Max power per module (W)' newline
                newline
                '%% Calculate PV power output' newline
                'P_dc = zeros(size(irradiance));' newline
                'for i = 1:length(irradiance)' newline
                '    % Irradiance effect' newline
                '    I_irr = I_sc * irradiance(i) / STC_irradiance;' newline
                '    ' newline
                '    % Temperature effect on Voc' newline
                '    T_cell = temperature(i) + (irradiance(i)/1000) * 25;  %% Nominal Operating Cell Temp' newline
                '    deltaT = T_cell - T_ref;' newline
                '    V_oc_T = V_oc * (1 + tempCoeff * deltaT);' newline
                '    ' newline
                '    % Simplified power calculation' newline
                '    P_per_module = P_mp * (irradiance(i) / STC_irradiance) * (1 + tempCoeff * deltaT);' newline
                '    P_dc(i) = numModules * P_per_module;' newline
                'end' newline
                newline
                '%% AC Power Output (with inverter loss)' newline
                'inverterEfficiency = 0.97;  %% Modern inverter efficiency' newline
                'P_ac = P_dc * inverterEfficiency / 1000;  %% Convert to MW' newline
                newline
                '%% Plot results' newline
                'figure;' newline
                'subplot(2,2,1); plot(hours, irradiance); xlabel(''Hour''); ylabel(''Irradiance (W/m²)''); title(''Solar Irradiance Profile'');' newline
                'subplot(2,2,2); plot(hours, temperature); xlabel(''Hour''); ylabel(''Temperature (°C)''); title(''Temperature Variation'');' newline
                'subplot(2,2,3); plot(hours, P_dc/1000); xlabel(''Hour''); ylabel(''DC Power (MW)''); title(''PV Array DC Output'');' newline
                'subplot(2,2,4); plot(hours, P_ac); xlabel(''Hour''); ylabel(''AC Power (MW)''); title(''Inverter AC Output'');' newline
                'sgtitle(''Solar PV System Performance'');' newline
                newline
                '%% Key metrics' newline
                'fprintf(''\\n=== Solar PV System Performance ===\\n'');' newline
                'fprintf(''Peak Power: %.2f MW\\n'', max(P_ac));' newline
                'fprintf(''Average Power: %.2f MW\\n'', mean(P_ac));' newline
                'fprintf(''Daily Energy: %.2f MWh\\n'', trapz(hours, P_ac)/24);' newline
                'fprintf(''Capacity Factor: %.2f %%\\n'', 100*mean(P_ac)/capacity);' newline
                newline
            ];
        end

        function code = getWindCode(obj)
            % Generate wind farm specific code
            code = [
                '%% Wind Farm Model' newline
                newline
                'numTurbines = %d;' newline, obj.weatherData.numTurbines
                'ratedPower = %.2f;  %% MW per turbine' newline, obj.weatherData.ratedPower
                'hubHeight = %.0f;  %% meters' newline, obj.weatherData.hubHeight
                'rotorDiameter = %.0f;  %% meters' newline, obj.weatherData.rotorDiameter
                'v_min = %.1f;  %% m/s - cut-in speed' newline, obj.weatherData.minWindSpeed
                'v_max = %.1f;  %% m/s - cut-out speed' newline, obj.weatherData.maxWindSpeed
                'v_r = 12;  %% m/s - rated wind speed' newline
                newline
                '%% Turbine Power Curve' newline
                '% Using a simplified power curve model' newline
                'v_range = 0:25;  %% Wind speed range (m/s)' newline
                'P_turbine = zeros(size(v_range));' newline
                'Cp_max = 0.45;  %% Max power coefficient' newline
                'A_rotor = pi * (rotorDiameter/2)^2;  %% Rotor area' newline
                'rho = 1.225;  %% Air density (kg/m³)' newline
                newline
                'for i = 1:length(v_range)' newline
                '    v = v_range(i);' newline
                '    if v < v_min || v > v_max' newline
                '        P_turbine(i) = 0;' newline
                '    elseif v <= v_r' newline
                '        % Cubic region (ramp up)' newline
                '        P_turbine(i) = ratedPower * (v - v_min)^3 / (v_r - v_min)^3;' newline
                '    else' newline
                '        % Constant power region' newline
                '        P_turbine(i) = ratedPower;' newline
                '    end' newline
                'end' newline
                newline
                '%% Wind Speed Time Series (Example)' newline
                '% Rayleigh distribution wind speed' newline
                 'time = 0:0.1:24;  %% 24 hours' newline
                'v_mean = 9;  %% Mean wind speed (m/s)' newline
                'v_wind = max(0, v_mean * sqrt(-2*log(rand(size(time)))));  %% Rayleigh samples' newline
                '% Apply low-pass filter for realistic variation' newline
                'v_wind = filtfilt(ones(1,10)/10, 1, v_wind);' newline
                newline
                '%% Wind Farm Power Output' newline
                'P_farm = zeros(size(v_wind));' newline
                'wake_loss = 0.1;  %% 10% wake losses' newline
                'for i = 1:length(v_wind)' newline
                '    % Interpolate power curve' newline
                '    P_per_turbine = interp1(v_range, P_turbine, v_wind(i), ''linear'', 0);' newline
                '    P_farm(i) = numTurbines * P_per_turbine * (1 - wake_loss);' newline
                'end' newline
                'P_farm = P_farm / 1000;  %% Convert to MW' newline
                newline
                '%% Plot results' newline
                'figure;' newline
                'subplot(2,2,1); plot(v_range, P_turbine); xlabel(''Wind Speed (m/s)''); ylabel(''Power (MW)''); title(''Turbine Power Curve''); grid on;' newline
                'subplot(2,2,2); plot(time, v_wind); xlabel(''Time (hours)''); ylabel(''Wind Speed (m/s)''); title(''Wind Speed Profile''); grid on;' newline
                'subplot(2,2,3); plot(time, P_farm); xlabel(''Time (hours)''); ylabel(''Farm Power (MW)''); title(''Wind Farm AC Output''); grid on;' newline
                'subplot(2,2,4); histogram(v_wind, 20); xlabel(''Wind Speed (m/s)''); ylabel(''Frequency''); title(''Wind Speed Distribution'');' newline
                'sgtitle(''Wind Farm Analysis'');' newline
                newline
                '%% Key metrics' newline
                'fprintf(''\\n=== Wind Farm Performance ===\\n'');' newline
                'fprintf(''Peak Power: %.2f MW\\n'', max(P_farm));' newline
                'fprintf(''Average Power: %.2f MW\\n'', mean(P_farm));' newline
                'fprintf(''Energy (24h): %.2f MWh\\n'', trapz(time, P_farm));' newline
                'fprintf(''Capacity Factor: %.2f %%\\n'', 100*mean(P_farm)/(numTurbines*ratedPower));' newline
                'fprintf(''Mean Wind Speed: %.2f m/s\\n'', mean(v_wind));' newline
                newline
            ];
        end

        function code = getStorageCode(obj)
            % Generate battery storage specific code
            code = [
                '%% Battery Energy Storage System (BESS)' newline
                newline
                'capacity_MWh = %.2f;  %% Energy capacity' newline, obj.storageCapacity
                'power_MW = %.2f;  %% Power rating' newline, obj.capacity
                'eta_round_trip = 0.85;  %% Round-trip efficiency' newline
                'charge_time = 2;  %% hours' newline
                'discharge_time = 2;  %% hours' newline
                'init_SOC = 50;  %% Initial state of charge (%)' newline
                newline
                '%% Battery Charging/Discharging Profile' newline
                'time = 0:0.1:24;  %% 24 hours' newline
                'P_demand = zeros(size(time));' newline
                'P_supply = zeros(size(time));' newline
                newline
                '%% Create demand profile' newline
                'for t = 1:length(time)' newline
                '    hour = mod(time(t), 24);' newline
                '    % Peak demand: 8-10am, 6-8pm' newline
                '    if (hour >= 8 && hour <= 10) || (hour >= 18 && hour <= 20)' newline
                '        P_demand(t) = 8;  %% MW' newline
                '    elseif hour >= 13 && hour <= 15' newline
                '        P_demand(t) = 5;  %% MW (shoulder)' newline
                '    else' newline
                '        P_demand(t) = 2;  %% MW (baseload)' newline
                '    end' newline
                'end' newline
                newline
                '%% Battery dispatch with simple rule-based control' newline
                'SOC = init_SOC * ones(size(time));' newline
                'P_batt = zeros(size(time));  %% Battery power (+ = discharge, - = charge)' newline
                'dt = time(2) - time(1);  %% timestep in hours' newline
                newline
                'for t = 2:length(time)' newline
                '    hour = mod(time(t), 24);' newline
                '    ' newline
                '    % Simple dispatch: charge during off-peak, discharge at peak' newline
                '    if (hour >= 0 && hour < 6) && SOC(t-1) < 90' newline
                '        % Charging window - charge from grid' newline
                '        P_batt(t) = -power_MW;  %% Charging power' newline
                '    elseif (hour >= 18 && hour <= 22) && SOC(t-1) > 20' newline
                '        % Discharging window - supply peak demand' newline
                '        P_batt(t) = power_MW;  %% Discharging power' newline
                '    else' newline
                '        P_batt(t) = 0;' newline
                '    end' newline
                '    ' newline
                '    % Update SOC' newline
                '    SOC(t) = SOC(t-1);' newline
                '    if P_batt(t) < 0  %% Charging' newline
                '        energy_change = -P_batt(t) * dt * eta_round_trip^0.5;' newline
                '    else  %% Discharging' newline
                '        energy_change = P_batt(t) * dt / eta_round_trip^0.5;' newline
                '    end' newline
                '    SOC(t) = SOC(t) - energy_change * 100 / (capacity_MWh * power_MW);' newline
                '    ' newline
                '    % Limit SOC to 10-90%% operational range' newline
                '    SOC(t) = max(10, min(90, SOC(t)));' newline
                'end' newline
                newline
                '%% Calculate net power (after battery dispatch)' newline
                'P_net = P_demand + P_batt;  %% Net power needed from grid' newline
                newline
                '%% Plot results' newline
                'figure;' newline
                'subplot(3,1,1); plot(time, P_demand, ''k-'', time, P_batt, ''b-''); xlabel(''Time (hours)''); ylabel(''Power (MW)''); title(''Demand and Battery Power''); legend(''Demand'', ''Battery'');' newline
                'subplot(3,1,2); plot(time, P_net, ''r-''); xlabel(''Time (hours)''); ylabel(''Grid Power (MW)''); title(''Net Power to Grid'');' newline
                'subplot(3,1,3); plot(time, SOC, ''g-''); xlabel(''Time (hours)''); ylabel(''SOC (%%)''); title(''Battery State of Charge''); ylim([0 100]);' newline
                'sgtitle(''Battery Energy Storage System Analysis'');' newline
                newline
                '%% Key metrics' newline
                'fprintf(''\\n=== BESS Performance ===\\n'');' newline
                'fprintf(''Peak Discharge Power: %.2f MW\\n'', max(P_batt));' newline
                'fprintf(''Energy Dispatched: %.2f MWh\\n'', sum(P_batt(P_batt>0)) * (time(2)-time(1)));' newline
                'fprintf(''Peak Reduction: %.2f %%\\n'', 100*(max(P_demand)-max(P_net))/max(P_demand));' newline
                'fprintf(''Final SOC: %.1f %%\\n'', SOC(end));' newline
                newline
            ];
        end

        function code = getHybridCode(obj)
            % Generate hybrid PV+Wind system code
            code = [
                '%% Hybrid PV + Wind System' newline
                newline
                'solar_capacity = %.2f;  %% MW' newline, obj.weatherData.solarCapacity
                'wind_capacity = %.2f;  %% MW' newline, obj.weatherData.windCapacity
                'storage_capacity = %.2f;  %% MWh' newline, obj.weatherData.storageCapacity
                newline
                '%% Combined renewable generation profile' newline
                'time = 0:0.1:24;' newline
                newline
                '%% Solar generation (peak at noon)' newline
                'P_solar = solar_capacity * max(0, sin(pi * time / 24)).^1.5;' newline
                newline
                '%% Wind generation (complementary to solar)' newline
                'v_wind = 7 + 4*sin(time/24 + pi/3);  %% Wind speed pattern' newline
                'P_wind = wind_capacity * (v_wind/12).^3;  %% Simplified power curve' newline
                'P_wind = min(wind_capacity, max(0, P_wind));' newline
                newline
                '%% Total renewable generation' newline
                'P_renewable = P_solar + P_wind;' newline
                newline
                '%% Storage dispatch (simple peak shaving)' newline
                'P_demand = 15 + 5*sin(pi*time/12);  %% Example demand' newline
                'P_stored = zeros(size(time));' newline
                'SOC = 50 * ones(size(time));' newline
                'dt = time(2) - time(1);' newline
                newline
                'for t = 2:length(time)' newline
                '    surplus = P_renewable(t) - P_demand(t);' newline
                '    if surplus > 0 && SOC(t-1) < 90' newline
                '        P_stored(t) = min(surplus, 50);  %% Charge' newline
                '    elseif surplus < 0 && SOC(t-1) > 20' newline
                '        P_stored(t) = max(surplus, -50);  %% Discharge' newline
                '    end' newline
                '    SOC(t) = SOC(t-1) - P_stored(t) * dt / storage_capacity;' newline
                '    SOC(t) = max(20, min(90, SOC(t)));' newline
                'end' newline
                newline
                '%% Grid power requirement' newline
                'P_grid = P_demand - P_renewable - P_stored;' newline
                newline
                '%% Plots' newline
                'figure;' newline
                'subplot(2,1,1); plot(time, P_solar, ''y-'', time, P_wind, ''b-'', time, P_renewable, ''k-'', time, P_demand, ''r--''); xlabel(''Time (hours)''); ylabel(''Power (MW)''); title(''Hybrid System Generation vs Demand''); legend(''Solar'', ''Wind'', ''Total Renewable'', ''Demand'');' newline
                'subplot(2,1,2); plot(time, P_grid, ''g-'', time, SOC, ''b--''); ylabel(''Power (MW)''); yyaxis right; ylabel(''SOC (%%)''); title(''Grid Power and Storage State'');' newline
                newline
                '%% Metrics' newline
                'fprintf(''Renewable Penetration: %.1f %%\\n'', 100*trapz(time, P_renewable)/(trapz(time, P_demand)+1e-6));' newline
                'fprintf(''Grid Import: %.1f MWh\\n'', trapz(time, max(0, P_grid)));' newline
                'fprintf(''Grid Export: %.1f MWh\\n'', trapz(time, max(0, -P_grid)));' newline
                newline
            ];
        end

        function code = getMicrogridCode(obj)
            % Generate microgrid system code
            code = [
                '%% Microgrid System Analysis' newline
                newline
                'numDGs = %d;  %% Number of DG units' newline, obj.weatherData.numDGs
                'totalLoad = %.2f;  %% MW' newline, obj.weatherData.totalLoad
                'controlMode = ''%s'';' newline, obj.weatherData.controlArchitecture
                'canIsland = %d;  %% Island capability' newline, obj.weatherData.canIsland
                newline
                '%% Distributed Generation (DG) Units' newline
                'DG = struct();' newline
                'for i = 1:numDGs' newline
                '    DG(i).capacity = totalLoad/numDGs;  %% MW' newline
                '    DG(i).type = mod(i, 3) + 1;  %% Types: 1=Diesel, 2=Renewable, 3=Storage' newline
                '    DG(i).voltage = 0.48;  %% kV' newline
                '    DG(i).status = 1;  %% Online' newline
                '    DG(i).Pgen = 0;' newline
                '    DG(i).Qgen = 0;' newline
                'end' newline
                newline
                '%% Demand Profile' newline
                'time = 0:0.25:24;  %% 15-minute intervals' newline
                'hours = mod(time, 24);' newline
                'load = totalLoad * (0.5 + 0.3*sin(pi*hours/12) + 0.1*sin(2*pi*hours/24));' newline
                'load = max(0.2*totalLoad, load);  %% Minimum 20%% load' newline
                newline
                '%% Grid-Connected Operation' newline
                'fprintf(''\\n=== Microgrid Grid-Connected Mode ===\\n'');' newline
                'P_grid = zeros(size(time));' newline
                'P_dg = zeros(numDGs, length(time));' newline
                'Q_total = zeros(size(time));' newline
                newline
                'for t = 1:length(time)' newline
                '    % Central dispatch (centralized control)' newline
                '    available_gen = sum([DG.capacity]);' newline
                '    ' newline
                '    if strcmp(controlMode, ''centralized'')' newline
                '        % Dispatch DG proportionally' newline
                '        for i = 1:numDGs' newline
                '            P_dg(i,t) = load(t) * DG(i).capacity / available_gen;' newline
                '        end' newline
                '    else' newline
                '        % Decentralized droop control' newline
                '        for i = 1:numDGs' newline
                '            P_dg(i,t) = DG(i).capacity * max(0, 1 - (1.05 - 0.95) / (2*DG(i).capacity));' newline
                '        end' newline
                '    end' newline
                '    ' newline
                '    % Reactive power control' newline
                '    Q_total(t) = load(t) * 0.3;  %% Assume 0.3 pu reactive load' newline
                '    ' newline
                '    % Power balance' newline
                '    P_grid(t) = load(t) - sum(P_dg(:,t));' newline
                'end' newline
                newline
                '%% Island Mode Analysis' newline
                'if canIsland' newline
                '    fprintf(''\\n=== Microgrid Island Mode ===\\n'');' newline
                '    ' newline
                '    % Frequency response during transient' newline
                '    t_island = 5;  %% Island command at t=5s' newline
                '    freq_nominal = 60;  %% Hz' newline
                '    freq = freq_nominal * ones(size(time));' newline
                '    m = 2;  %% Inertia coefficient' newline
                '    d = 0.1;  %% Damping' newline
                '    ' newline
                '    for t = 2:length(time)' newline
                '        if time(t) >= t_island' newline
                '            % Frequency deviation based on power balance' newline
                '            dP = load(t) - sum(P_dg(:,t));  % Imbalance' newline
                '            df_dt = -dP / (m * totalLoad) - d * (freq(t-1) - freq_nominal);' newline
                '            freq(t) = freq(t-1) + df_dt * (time(2)-time(1));' newline
                '        end' newline
                '    end' newline
                'end' newline
                newline
                '%% Results and Plotting' newline
                'figure;' newline
                'subplot(2,2,1); plot(time, load, ''k-'', ''LineWidth'', 2); hold on;' newline
                'for i = 1:numDGs; plot(time, P_dg(i,:)); end; xlabel(''Time (hours)''); ylabel(''Power (MW)''); title(''DG Dispatch''); legend([{''Load''}, cell(numDGs,1)]);' newline
                'subplot(2,2,2); plot(time, P_grid); xlabel(''Time (hours)''); ylabel(''Grid Power (MW)''); title(''Grid Power Exchange'');' newline
                'if canIsland; subplot(2,2,3); plot(time, freq); xlabel(''Time (hours)''); ylabel(''Frequency (Hz)''); title(''Frequency Response''); yline(60, ''k--''); end' newline
                'subplot(2,2,4); plot(time, Q_total); xlabel(''Time (hours)''); ylabel(''Reactive Power (MVAr)''); title(''Reactive Power'');' newline
                'sgtitle(''Microgrid Analysis'');' newline
                newline
                '%% Summary' newline
                'fprintf(''Peak Load: %.2f MW\\n'', max(load));' newline
                'fprintf(''Total Energy Demand: %.2f MWh\\n'', trapz(time, load));' newline
                'fprintf(''Grid Import: %.2f MWh\\n'', trapz(time, max(0, P_grid)));' newline
                'fprintf(''Grid Export: %.2f MWh\\n'', trapz(time, max(0, -P_grid)));' newline
                newline
            ];
        end

        function code = getCommonAnalysisCode(obj)
            % Generate common analysis code for all renewable types
            code = [
                '%% Control Performance Analysis' newline
                'fprintf(''\\n=== Control Mode: %s ===\\n'', controlMode);' newline
                newline
                'switch controlMode' newline
                '    case ''PQ''' newline
                '        fprintf(''Constant Power Control: P and Q are independently specified\\n'');' newline
                '        fprintf(''Suitable for: Grid-following converters\\n'');' newline
                '    case ''PV''' newline
                '        fprintf(''PV Control: Voltage magnitude is maintained at setpoint\\n'');' newline
                '        fprintf(''Suitable for: Weak grid, voltage support\\n'');' newline
                '    case ''droop''' newline
                '        fprintf(''Droop Control: Frequency and voltage droop characteristics\\n'');' newline
                '        fprintf(''Suitable for: Grid-supporting operation\\n'');' newline
                '    case ''grid_forming''' newline
                '        fprintf(''Grid-Forming Control: Full voltage control capability\\n'');' newline
                '        fprintf(''Suitable for: Strong grid support, fault ride-through\\n'');' newline
                'end' newline
                newline
                '%% Penetration Level Impact' newline
                'fprintf(''\\nPenetration Level: %.1f %%\\n'', penetrationLevel);' newline
                'if penetrationLevel > 50' newline
                '    fprintf(''** High penetration level - advanced control required **\\n'');' newline
                '    fprintf(''Recommended: Grid-forming converters with FIRT capability\\n'');' newline
                'elseif penetrationLevel > 20' newline
                '    fprintf(''** Medium penetration - coordinated control needed **\\n'');' newline
                '    fprintf(''Recommended: Grid-supporting controls with ramp rate limits\\n'');' newline
                'else' newline
                '    fprintf(''** Low penetration - standard grid-following acceptable **\\n'');' newline
                'end' newline
                newline
                '%% Harmonic Vulnerability Assessment' newline
                'fprintf(''\\n=== Harmonic Assessment ===\\n'');' newline
                'short_circuit_ratio = 100 / (penetrationLevel + 1);  %% Simplified SCR' newline
                'fprintf(''Short Circuit Ratio: %.2f\\n'', short_circuit_ratio);' newline
                'if short_circuit_ratio < 3' newline
                '    fprintf(''Warning: Weak grid conditions - elevated harmonic distortion expected\\n'');' newline
                'end' newline
                newline
            ];
        end

        function handleInput(obj, userInput)
            % Handle contextual user input
            fprintf('Renewable integration module received: %s\n', userInput);
        end
    end
end
