%% Stability Studies - Transient Stability Analysis Template
% Auto-generated MATLAB code for power system transient stability analysis
%
% This script performs time-domain simulation of power system dynamics
% following a disturbance (fault, line trip, or generator trip).

clear; clc; close all;

%% System Configuration
baseMVA = 100;
systemFrequency = 60;  % Hz
simulationTime = 15;   % seconds
timeStep = 0.001;      % 1 ms time step

% ODE Solver selection
solverType = 'ode23t';  % Good for stiff power system equations

fprintf('\n================================================================================\n');
fprintf('           Power System Transient Stability Analysis\n');
fprintf('================================================================================\n\n');

%% Generator Parameters (Example: 2-machine system)
gen(1).Sn = 900;       % Rated power (MVA)
gen(1).Vn = 20;        % Rated voltage (kV)
gen(1).H = 5.0;        % Inertia constant (s)
gen(1).D = 2.0;        % Damping coefficient
gen(1).wp = 1.0;       % Initial operating point (pu)

gen(2).Sn = 900;
gen(2).Vn = 20;
gen(2).H = 4.5;
gen(2).D = 1.8;
gen(2).wp = 1.0;

%% Network Data
% 2-machine test system
Zline = 0.05 + 1i*0.25;  % Line impedance (pu)

%% Disturbance Definition
disturbanceType = 'three_phase_fault';
disturbanceTime = 1.0;   % seconds
faultDuration = 0.1;     % 100 ms fault clearing time
faultLocation = 0.5;     % Middle of line between bus 1 and 2

fprintf('Disturbance: %s at t = %.2f s\n', disturbanceType, disturbanceTime);
fprintf('Fault clearing: t = %.2f s (duration: %.3f s)\n\n', ...
    disturbanceTime + faultDuration, faultDuration);

%% Initial Conditions
n_gen = length(gen);
delta = zeros(n_gen, 1);    % Rotor angles (rad)
omega = ones(n_gen, 1);     % Angular velocities (pu)
Pm = ones(n_gen, 1) * 0.9;  % Mechanical power (pu)
Eq = ones(n_gen, 1);        % Transient EMF (pu)

% State vector: [delta1, omega1, delta2, omega2, ...]
x0 = [delta; omega];

%% ODE Solver
options = odeset('MaxStep', timeStep, 'RelTol', 1e-6, 'AbsTol', 1e-8);
tspan = [0, simulationTime];

fprintf('Solving differential equations (this may take a moment)...\n');

% Use ode23t for stiff system
[t, x] = ode23t(@(t,x) powerSystemDynamics(t, x, gen, baseMVA, ...
    disturbanceTime, faultDuration, Pm, Zline), tspan, x0, options);

fprintf('Solution complete. Simulated %d time steps.\n\n', length(t));

%% Extract results
delta_1 = x(:, 1);
omega_1 = x(:, 2);
if n_gen > 1
    delta_2 = x(:, 3);
    omega_2 = x(:, 4);
end

% Rotor angle difference
delta_diff = delta_1 - delta_2;

%% Stability Assessment
omega_min = min(omega_1);
omega_max = max(omega_1);
delta_max = max(abs(delta_diff));

fprintf('=== Stability Assessment ===\n');
fprintf('Generator 1:\n');
fprintf('  Min frequency: %.4f pu (%.2f Hz)\n', omega_min, omega_min*systemFrequency);
fprintf('  Max frequency: %.4f pu (%.2f Hz)\n', omega_max, omega_max*systemFrequency);
fprintf('  Frequency nadir: %.2f Hz at t=%.3f s\n', min(omega_1)*systemFrequency, t(find(omega_1 == min(omega_1), 1)));

if delta_diff(end) < 1.5  % 85 degrees
    fprintf('RESULT: System is STABLE\n');
    fprintf('Final rotor angle difference: %.2f degrees\n', delta_diff(end)*180/pi);
else
    fprintf('RESULT: System is UNSTABLE - generators lost synchronism\n');
    fprintf('Max rotor angle difference: %.2f degrees\n', delta_max*180/pi);
end
fprintf('\n');

%% Plots
figure('Position', [100, 100, 1200, 800]);

% Rotor angles
subplot(2,3,1);
plot(t, delta_1*180/pi, 'b-', 'LineWidth', 1.5);
if n_gen > 1
    hold on;
    plot(t, delta_2*180/pi, 'r-', 'LineWidth', 1.5);
    legend('Generator 1', 'Generator 2');
end
xlabel('Time (s)');
ylabel('Rotor Angle (degrees)');
title('Rotor Angle');
grid on;
xline(disturbanceTime, 'k--', 'Fault');
xline(disturbanceTime + faultDuration, 'k--', 'Cleared');

% Rotor angle difference
subplot(2,3,2);
plot(t, delta_diff*180/pi, 'g-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Angle Difference (degrees)');
title('Relative Rotor Angle (Gen1 - Gen2)');
grid on;
yline(0, 'k--');
yline(90, 'r--', 'Stability Limit');
xline(disturbanceTime, 'k--');
xline(disturbanceTime + faultDuration, 'k--');

% Frequency
subplot(2,3,3);
plot(t, omega_1*systemFrequency, 'b-', 'LineWidth', 1.5);
if n_gen > 1
    hold on;
    plot(t, omega_2*systemFrequency, 'r-', 'LineWidth', 1.5);
    legend('Generator 1', 'Generator 2');
end
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Generator Frequency');
grid on;
yline(systemFrequency, 'k--', 'Nominal');
xline(disturbanceTime, 'k--');
xline(disturbanceTime + faultDuration, 'k--');

% Angular velocity (pu)
subplot(2,3,4);
plot(t, omega_1, 'b-', 'LineWidth', 1.5);
if n_gen > 1
    hold on;
    plot(t, omega_2, 'r-', 'LineWidth', 1.5);
end
xlabel('Time (s)');
ylabel('Angular Velocity (pu)');
title('Generator Angular Velocity');
grid on;
yline(1.0, 'k--');
xline(disturbanceTime, 'k--');
xline(disturbanceTime + faultDuration, 'k--');

% Rate of change of frequency (RoCoF)
subplot(2,3,5);
rocof = gradient(omega_1, t) * systemFrequency;
plot(t, rocof, 'b-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('RoCoF (Hz/s)');
title('Rate of Change of Frequency');
grid on;
xline(disturbanceTime, 'k--');
xline(disturbanceTime + faultDuration, 'k--');

% Phase plane (delta vs omega)
subplot(2,3,6);
plot(delta_1*180/pi, omega_1, 'b-', 'LineWidth', 1.5);
if n_gen > 1
    hold on;
    plot(delta_2*180/pi, omega_2, 'r-', 'LineWidth', 1.5);
end
xlabel('Rotor Angle (degrees)');
ylabel('Angular Velocity (pu)');
title('Phase Plane - Rotor Dynamics');
grid on;
yline(1.0, 'k--');
legend('Generator 1', 'Generator 2');

sgtitle('Transient Stability Analysis Results', 'FontSize', 12, 'FontWeight', 'bold');

%% Power Flow Calculation (simplified)
fprintf('\nLine Power Flow (Average):\n');
fprintf('  Bus 1 to Bus 2: %.2f MW\n', mean(delta_1 - delta_2) / norm(Zline) * baseMVA);
fprintf('\n');

%% Function: Power System Dynamics
function dxdt = powerSystemDynamics(t, x, gen, baseMVA, ...
    disturbanceTime, faultDuration, Pm, Zline)
    
    n_gen = length(gen);
    
    % Extract states
    delta = x(1:n_gen);
    omega = x(n_gen+1:end);
    
    % Initialize derivatives
    dxdt = zeros(2*n_gen, 1);
    
    % Rotor angle rates
    dxdt(1:n_gen) = omega;
    
    % Electrical power (simplified)
    % For 2-machine system: Pe = (V1*V2/X) * sin(delta1 - delta2)
    if n_gen == 2
        X = imag(Zline);
        Pe = sin(delta(1) - delta(2)) / X;
        
        % Determine Pe during fault
        if t >= disturbanceTime && t < (disturbanceTime + faultDuration)
            Pe = Pe * 0.1;  % Reduced power transfer during fault
        end
    else
        Pe = 0;
    end
    
    % Swing equation: d(omega)/dt = (1/(2*H)) * (Pm - Pe - D*(omega-1))
    for i = 1:n_gen
        H = gen(i).H;
        D = gen(i).D;
        
        % Mechanical input power
        Pm_i = Pm(i);
        
        % Electrical output power (simplified)
        if i == 1
            Pe_i = Pe;
        else
            Pe_i = -Pe;
        end
        
        % Swing equation
        dxdt(n_gen + i) = (Pm_i - Pe_i - D*(omega(i) - 1)) / (2*H);
    end
end
