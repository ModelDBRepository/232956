%% Demonstration script for new NVU model

%% Version control
% This is the latest version of the NVU model: Version 1.2!
% K+, NO, Astrocytic Ca2+, TRPV4, ECS

%% Construct NVU
% The NVU consists of a number of submodules, implemented as MATLAB
% classes, presently an astrocyte, a lumped SMC/EC model, and a model of
% the wall mechanics. 
%
% The parameters of each of these submodules are
% specified when the modules are constructed, here in the call to NVU, see
% for example the |SMCEC| part of the |NVU| call below:
%
% Options for the ODE solver (currently |ode15s|) are provided by
% specifying the |odeopts| parameter. The code works fine with default
% tolerances.
clear all

odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 0.5, 'Vectorized', 1);

XLIM1 = 0; XLIM2 = 500;
FIG_NUM = 1;


NEURONAL_START  = 100;      % Start of neuronal stimulation
NEURONAL_END    = 300;      % End of neuronal stimulation 
ECS_START       = 100000000;      % Start of ECS K+ input
ECS_END         = 200000000;      % End of ECS K+ input
J_PLC           = 0.11;      % Jplc value in EC: 0.11 for steady state, 0.3 for oscillations

ECS             = 1;        % Include ECS compartment or not
CA_SWITCH       = 1;        % Turn on Ca2+ pathway
NO_INPUT_SWITCH = 1;        % Turn on NO stimulation
NO_PROD_SWITCH  = 1;        % Turn on Nitric Oxide production 
K_SWITCH        = 1;        % Turn on K+ input into SC
TRPV_SWITCH     = 1;        % Turn on TRPV4 Ca2+ channel from AC to PVS

nv = NVU(Neuron('F_input', 2.67, 'startpulse', NEURONAL_START, 'lengthpulse', NEURONAL_END - NEURONAL_START, 'GluSwitch', NO_INPUT_SWITCH, 'NOswitch', NO_PROD_SWITCH, 'KSwitch', K_SWITCH), ...
    Astrocyte('J_max', 2880, 'ECS_input', 9, 'trpv_switch', TRPV_SWITCH, 'startpulse', NEURONAL_START, 'lengthpulse', NEURONAL_END - NEURONAL_START, 't0_ECS', ECS_START, 'tend_ECS', ECS_END, 'GluSwitch', CA_SWITCH, 'ECSswitch', ECS, 'PVStoECS', 0, 'SCtoECS', 1), ...
    WallMechanics(), ...
    SMCEC('J_PLC', J_PLC, 'NOswitch', NO_PROD_SWITCH), 'odeopts', odeopts);

nv.T = linspace(0, XLIM2*2, 5000);
nv.simulate()

% Run this to take the ICs as the end of the last simulation run i.e. steady state ICs
ICs = (nv.U(end, :))';
nv.u0 = ICs;
nv.simulateManualICs() 

% Plot figures - whatever you want
figure(FIG_NUM);

subplot(2,2,1);
    hold all;
    plot(nv.T, nv.out('Ca_k'), 'LineWidth', 1);
    ylabel('Ca_k [\muM]');
    xlim([XLIM1 XLIM2])
subplot(2,2,2);
    hold all;
    plot(nv.T, nv.out('w_k'), 'LineWidth', 1);
    ylabel('w_k [-]');
    xlim([XLIM1 XLIM2])
subplot(2,2,3);
    hold all;
    plot(nv.T, nv.out('K_p')/1e3, 'LineWidth', 1);
    ylabel('K_p [mM]');
    xlim([XLIM1 XLIM2])
subplot(2,2,4);
    hold all;
    plot(nv.T, nv.out('v_k')*1e3, nv.T, nv.out('E_BK_k')*1e3, 'LineWidth', 1);
    ylabel('v_k and E_{BK}');
    xlim([XLIM1 XLIM2])

figure(FIG_NUM + 1);

subplot(2,3,1)
    hold all;
    plot(nv.T, nv.out('Ca_k'), 'LineWidth', 1);
    ylabel('Ca_k [\muM]');
    xlim([XLIM1 XLIM2])
    xlabel('Time [s]'); 
subplot(2,3,2)
    hold all;
    plot(nv.T, nv.out('eet_k'), 'LineWidth', 1);
    ylabel('eet_k [\muM]');
    xlim([XLIM1 XLIM2])
    xlabel('Time [s]'); 
subplot(2,3,3)
    hold all;
    plot(nv.T, nv.out('w_k'), 'LineWidth', 1);
    ylabel('w_k [-]');
    xlim([XLIM1 XLIM2])
    xlabel('Time [s]'); 
subplot(2,3,4)
    hold all;
    plot(nv.T, nv.out('R')*1e6, 'LineWidth', 1);
    ylabel('Radius [\mum]');
    xlim([XLIM1 XLIM2])
    xlabel('Time [s]'); 
subplot(2,3,5)
    hold all;
    plot(nv.T, nv.out('K_p')/1e3, 'LineWidth', 1);
    xlabel('Time [s]'); 
    ylabel('K_p [mM]');
    xlim([XLIM1 XLIM2])
subplot(2,3,6)
    hold all;
    plot(nv.T, nv.out('R')*1e6, 'LineWidth', 1)
    xlabel('Time [s]'); 
    ylabel('Radius [\mum]');
    xlim([XLIM1 XLIM2])
