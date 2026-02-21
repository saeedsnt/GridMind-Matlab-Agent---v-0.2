%% MATLAB AI Agent - Quick Reference Card
%
% ============================================================================
%                          QUICK START
% ============================================================================
%
% Start the agent:
%   main
%
% Available Commands:
%   /power-flow       - Load flow analysis
%   /power-fault      - Fault analysis
%   /power-stability  - Stability studies
%   /power-protection - Protection relay design
%   /power-renewable  - Renewable integration
%   /help             - Show this help
%   /quit             - Exit
%
% ============================================================================
%                      POWER FLOW ANALYSIS
% ============================================================================
%
% Command: /power-flow
% Best for: Steady-state voltage and power flow studies
%
% You will enter:
%   • Number of buses (2-10000): e.g., 5
%   • Base MVA (1-10000): e.g., 100
%   • Solution method: Newton-Raphson (recommended)
%   • For each bus:
%       - Bus type (PQ=load, PV=generator, slack=reference)
%       - Voltage/Power specifications
%   • For each line:
%       - Resistance, reactance, susceptance
%       - Tap ratio (transformer)
%
% Output: MATLAB script with voltage profiles, power flows, and losses
% Save location: Current directory as {date}_power_flow_script.m
%
% ============================================================================
%                       FAULT ANALYSIS
% ============================================================================
%
% Command: /power-fault
% Best for: Short circuit studies and breaker/relay ratings
%
% You will enter:
%   • Fault type:
%       - 3-phase (symmetrical)
%       - Single line-to-ground (SLG)
%       - Line-to-line (LL)
%       - Double line-to-ground (DLG)
%   • Fault location (bus number)
%   • Fault impedance (pu)
%   • System X/R ratio
%   • Pre-fault voltage (pu)
%
% Output: Fault current calculations, time-domain waveforms
% Includes: Symmetrical RMS, DC component, breaker duty
%
% ============================================================================
%                    TRANSIENT STABILITY
% ============================================================================
%
% Command: /power-stability
% Best for: Rotor angle stability and frequency response
%
% Study Types:
%   1. Transient Stability (time-domain simulation)
%   2. Small Signal Stability (eigenvalue analysis)
%   3. Voltage Stability (P-V curves)
%   4. Frequency Stability (load-frequency dynamics)
%
% Parameters Required:
%   • Simulation time (typically 10-20 seconds)
%   • Time step (0.001 s recommended)
%   • ODE solver (ode23t for power systems)
%   • Generator data (H constant, damping, rating)
%   • Disturbance type (fault, line trip, generator trip)
%   • Disturbance timing and duration
%
% Output: Rotor angle trajectories, frequency response, stability margin
% Assessment: STABLE or UNSTABLE with critical clearing time
%
% ============================================================================
%                    RENEWABLE INTEGRATION
% ============================================================================
%
% Command: /power-renewable
% Best for: Solar, wind, battery storage, and microgrid analysis
%
% Resource Types:
%   1. Solar PV Systems
%   2. Wind Farms
%   3. Battery Energy Storage (BESS)
%   4. Hybrid PV + Wind
%   5. Microgrids
%
% For Solar:
%   • Number of modules/array size
%   • Temperature coefficient
%   • Module rating, efficiency
%   • Reference conditions
%
% For Wind:
%   • Number of turbines
%   • Turbine rated power
%   • Hub height, rotor diameter
%   • Min/max wind speeds
%   • Generator type (DFIG, PMSG, sync)
%
% For Storage:
%   • Energy capacity (MWh)
%   • Power rating (MW)
%   • Chemistry (Li-ion, lead-acid, flow)
%   • Round-trip efficiency
%   • Charge/discharge time
%
% Output: Generation profiles, grid impact, control performance
%
% ============================================================================
%                      PROTECTION SYSTEMS
% ============================================================================
%
% Command: /power-protection
% Best for: Relay coordination and protection scheme design
%
% Relay Functions Supported:
%   • 50/51 - Instantaneous/Time-delay overcurrent
%   • 21 - Distance (impedance) protection
%   • 87 - Differential protection
%   • 32 - Directional power protection
%
% You will specify:
%   • CT ratio (current transformer)
%   • VT ratio (voltage transformer)
%   • Pickup current (primary amperes)
%   • Time dial setting (TDS) [0.5-15]
%   • Curve type (IEEE, IEC standard)
%   • Coordination time interval (CTI)
%
% Output: Time-Current Characteristic (TCC) plots
% Shows: Coordination selectivity, relay sensitivity checks
%
% ============================================================================
%                    WORKFLOW EXAMPLE
% ============================================================================
%
% Session 1: Power Flow Study
%   >> main
%   > /power-flow
%   > 1 (select code-based)
%   > Enter 5 buses
%   > Enter 100 MVA base
%   > [Follow prompts]
%   > Save to file? y
%   > Filename: [default accepted]
%   > System generated power flow analysis script
%   > Run it: >> power_flow_script
%
% Session 2: Fault Analysis
%   > /power-fault
%   > 2 (select simulator)
%   > 1 (three-phase fault)
%   > [Continue workflow]
%
% ============================================================================
%                  KEYBOARD SHORTCUTS
% ============================================================================
%
% During parameter input:
%   ENTER        - Accept default value (if shown)
%   Ctrl+C       - Interrupt and return to command selection
%   'y' or 'Y'   - Confirm yes/no questions
%   'n' or 'N'   - Decline yes/no questions
%   [number]     - Select menu option by number
%
% ============================================================================
%                 PERFORMANCE TIPS
% ============================================================================
%
% For Fast Convergence:
%   • Start with fewer buses (5-10)
%   • Use flat start initialization
%   • Set tolerance = 1e-6 (tighter for accuracy)
%   • Use Newton-Raphson for best performance
%
% For Stability Analysis:
%   • Use time step = 0.001 s (1 ms)
%   • Use ode23t solver for stiff equations
%   • Set max simulation time = 15-20 seconds
%   • Use generator H = 3-8 seconds (typical)
%
% For Large Systems:
%   • Use Fast Decoupled Power Flow method
%   • Reduce time step to 0.01 s for faster simulation
%   • Use DC power flow for preliminary analysis
%   • Break systems into subsystems
%
% ============================================================================
%                    STANDARD VALUES (IEEE)
% ============================================================================
%
% Generator Inertia Constant (H):
%   Hydro Plant:        2-9 seconds
%   Thermal Plant:      3-6 seconds
%   Gas Turbine:        3-5 seconds
%   DFIG Wind:          2-4 seconds (synthetic)
%   Nuclear:            3-6 seconds
%
% System X/R Ratio:
%   Transmission:       10-20
%   Subtransmission:    5-15
%   Distribution:       1-5
%
% Voltage Limits (normal):
%   Transmission:       ±5% (0.95-1.05 pu)
%   Distribution:       ±10% (0.90-1.10 pu)
%
% Fault Impedance:
%   Solid ground:       0.00-0.01 pu
%   Arc fault:          0.01-0.1 pu
%   High-resistance:    0.1-1.0 pu
%
% ============================================================================
%                     TROUBLESHOOTING
% ============================================================================
%
% Problem: "Class not found"
% Solution: Ensure all .m files are in correct directories
%           Run: addpath(genpath('.'))
%
% Problem: Power flow doesn't converge
% Solution: 1. Reduce tolerance to 1e-4
%           2. Try Gauss-Seidel method
%           3. Check system data validity
%           4. Increase max iterations to 100
%
% Problem: Stability shows unrealistic results
% Solution: 1. Verify time step (should be 0.001 s)
%           2. Check generator H constant (should be 2-10 s)
%           3. Verify fault duration (0.05-0.15 s typical)
%           4. Use ode23t solver
%
% ============================================================================
%                       EXPORT OPTIONS
% ============================================================================
%
% All generated code can be:
%   • Run directly in MATLAB
%   • Saved as .m files for later use
%   • Modified and extended
%   • Integrated into larger projects
%   • Used with Simulink models
%   • Parameters exported to Excel (optional)
%
% ============================================================================
%
% For full documentation: help PROJECT_SUMMARY
% For detailed module info: help {ModuleName}
% Example: help PowerFlowModule
% Example: help StabilityModule
%
% ============================================================================
