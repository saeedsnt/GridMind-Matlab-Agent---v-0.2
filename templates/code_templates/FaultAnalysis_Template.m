%% Fault Analysis Template
% Auto-generated MATLAB code for short circuit and fault analysis
%
% Modify parameters below and run to perform fault calculations

clear; clc; close all;

%% System Parameters
baseMVA = 100;
systemVoltage = 138;  % kV
XR_ratio = 10;        % System X/R ratio
faultType = 'three_phase';
faultLocation = 1;    % Bus number
faultImpedance = 0.01; % Fault impedance (pu)
prefaultVoltage = 1.0; % Pre-fault voltage (pu)

fprintf('\n================================================================================\n');
fprintf('              Fault Analysis - %s Fault\n', strtrim(faultType));
fprintf('================================================================================\n\n');

%% Generate sequence impedances
% Based on pre-faultpower flow
Zpos = 0.10;  % Positive sequence impedance (pu)
Zneg = 0.08;  % Negative sequence impedance (pu)
Zzero = 0.05; % Zero sequence impedance (pu) (lower for solidly grounded)

% Fault impedance
Zf = faultImpedance * (1 + 1i*XR_ratio) / sqrt(1 + XR_ratio^2);

%% Calculate fault currents
switch lower(faultType)
    case 'three_phase'
        % Three-phase symmetrical fault
        I_f = prefaultVoltage / (Zpos + Zf);
        fprintf('Three-Phase Fault Current: %.2f kA (per-unit: %.4f)\n', ...
            abs(I_f)*baseMVA/sqrt(3)/systemVoltage, abs(I_f));
        fprintf('Fault impedance: %.4f pu\n', abs(Zf));
        
        % DC component
        tau_dc = (Zpos + Zf) / (2*pi*60*0.0001);  % Time constant (simplified)
        fprintf('Time constant: %.4f s\n', tau_dc);
        
    case 'single_line_ground'
        % Single line-to-ground fault
        I_0 = prefaultVoltage / (Zpos + Zneg + Zzero + 3*Zf);
        I_f = 3 * I_0;  % Fault current
        fprintf('Single-Line-to-Ground Fault Current: %.2f kA\n', ...
            abs(I_f)*baseMVA/sqrt(3)/systemVoltage);
        fprintf('Zero sequence component: %.4f pu\n', abs(I_0));
        
    case 'line_to_line'
        % Line-to-line fault
        I_f = prefaultVoltage / (Zpos + Zneg + Zf);
        fprintf('Line-to-Line Fault Current: %.2f kA\n', ...
            abs(I_f)*sqrt(3)*baseMVA/sqrt(3)/systemVoltage);
        fprintf('Phase shift: %.1f degrees\n', angle(I_f)*180/pi);
        
    case 'double_line_ground'
        % Double line-to-ground fault
        I_f = prefaultVoltage * (Zneg + Zzero) / ((Zneg + Zf) * (Zzero + Zf) / (Zneg + Zzero + 2*Zf) + Zf);
        fprintf('Double-Line-to-Ground Fault Current: %.2f kA\n', ...
            abs(I_f)*baseMVA/sqrt(3)/systemVoltage);
end

fprintf('\n');

%% Breaker duty calculation
fprintf('Breaker Duty Calculation:\n');
fprintf('  Symmetrical fault current: %.2f kA\n', abs(I_f)*100);
fprintf('  System base: %.0f MVA at %.0f kV\n', baseMVA, systemVoltage);
fprintf('  Breaker contact parting time: 2-8 cycles\n');
fprintf('\n');

%% Transient current calculation
fprintf('Transient Current Profile:\n');
time = linspace(0, 0.1, 1000);  % 100 ms
I_ac = zeros(size(time));
I_dc = zeros(size(time));

tau = 0.01;  % Transient time constant
for i = 1:length(time)
    I_ac(i) = abs(I_f);  % AC component
    I_dc(i) = abs(I_f) * exp(-time(i)/tau);  % DC component
end
I_total = I_ac + I_dc;

%% Plot results
figure;
subplot(2,1,1);
plot(time*1000, I_total, 'LineWidth', 2);
hold on;
plot(time*1000, I_ac, '--', 'LineWidth', 1);
plot(time*1000, I_dc, '-.', 'LineWidth', 1);
xlabel('Time (ms)');
ylabel('Fault Current (pu)');
title('Fault Current vs Time');
legend('Total', 'AC Component', 'DC Component');
grid on;

% RMS current
I_rms = sqrt(I_ac.^2 + 0.5*I_dc.^2);
subplot(2,1,2);
plot(time*1000, I_rms, 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('RMS Current (pu)');
title('RMS Fault Current');
grid on;

sgtitle('Fault Analysis Results');

%% Relay coordination guidance
fprintf('\nRelay Coordination Guidance:\n');
fprintf('  CT Ratio: %.0f (common choices: 100, 200, 300, 400, 500)\n', baseMVA/10);
fprintf('  Pickup current: %.1f A (2-4x normal load)\n', 5);
fprintf('  Recommended TDS: 0.5-2.0 (inverse curve)\n');
fprintf('  Coordination time interval (CTI): %.1f cycles\n', 3);
fprintf('\n');
