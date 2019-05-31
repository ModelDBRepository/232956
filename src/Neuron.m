classdef Neuron < handle
    properties
        params
        u0
        index
        n_out
        idx_out
        enabled
    end
    methods
        function self = Neuron(varargin)
            self.params = parse_inputs(varargin{:});
            self.index = indices();
            self.u0 = initial_conditions(self.index);
            self.enabled = true(size(self.u0));
            [self.idx_out, self.n_out] = output_indices();
        end
        function [du, varargout] = rhs(self, t, u, NO_k)
            t = t(:).';
            p = self.params;
            idx = self.index;
            
            % NO pathway 
            Ca_n = u(idx.Ca_n, :);                  
            nNOS_act_n = u(idx.nNOS_act_n, :);
            NO_n = u(idx.NO_n, :);
            
            % K+ input
            J_NaK_n = p.k_C * self.input_K(t);
            
            % NO pathway
            w_NR2A = self.input_Glu(t) ./ (p.K_mA + self.input_Glu(t)); %[-] 
            w_NR2B = self.input_Glu(t) ./ (p.K_mB + self.input_Glu(t)); %[-]
            
            I_Ca = (-4 * p.v_n * p.G_M * p.P_Ca_P_M * (p.Ca_ex / p.M)) / (1 + exp(-80 * (p.v_n + 0.02))) * (exp(2 * p.v_n * p.F / (p.R_gas * p.T))) / (1 - exp(2 * p.v_n * p.F / (p.R_gas * p.T))); %[fA]
            I_Ca_tot = I_Ca .* (p.n_NR2A * w_NR2A + p.n_NR2B * w_NR2B); %[fA]
            
            CaM = Ca_n / p.m_c; %[uM]
            tau_nk = p.x_nk ^ 2 ./  (2 * p.D_cNO); %[s]
            
            p_NO_n = p.NOswitch * ( nNOS_act_n * p.V_max_NO_n * p.O2_n / (p.K_mO2_n + p.O2_n) * p.LArg_n / (p.K_mArg_n + p.LArg_n) ); %[uM/s]
            c_NO_n = p.k_O2_n * NO_n.^2 * p.O2_n; %[uM/s]
            d_NO_n = (NO_k - NO_n) ./ tau_nk; %[uM/s]
            
            du(idx.Ca_n, :) = (I_Ca_tot / (2 * p.F * p.V_spine) - (p.k_ex * (Ca_n - p.Ca_rest))) / (1 + p.lambda_buf); %[uM/s]
            du(idx.nNOS_act_n, :) = p.V_maxNOS * CaM ./ (p.K_actNOS + CaM) - p.mu2_n * nNOS_act_n; %[uM/s]
            du(idx.NO_n, :) = p_NO_n - c_NO_n + d_NO_n; %[uM/s]
            
            
            du = bsxfun(@times, self.enabled, du);
            if nargout == 2
            	Uout = zeros(self.n_out, size(u, 2));
            	Uout(self.idx_out.ft, :) = self.input_K(t);
            	Uout(self.idx_out.J_Na_n, :) = J_NaK_n;
                Uout(self.idx_out.Glu, :) = self.input_Glu(t);
            	Uout(self.idx_out.w_NR2A, :) = w_NR2A;
            	Uout(self.idx_out.w_NR2B, :) = w_NR2B;
            	Uout(self.idx_out.I_Ca, :) = I_Ca;
            	Uout(self.idx_out.I_Ca_tot, :) = I_Ca_tot;
            	Uout(self.idx_out.CaM, :) = CaM;
            	Uout(self.idx_out.tau_nk, :) = tau_nk;
            	Uout(self.idx_out.p_NO_n, :) = p_NO_n;
            	Uout(self.idx_out.c_NO_n, :) = c_NO_n;
            	Uout(self.idx_out.d_NO_n, :) = d_NO_n;
           
            	varargout = {Uout};
            end
        end
        
       function [J_Na_n, NO_n] = shared(self, t, u)    %shared variables
            t = t(:).';
            p = self.params;
            idx = self.index;
            J_Na_n = p.k_C * self.input_K(t);
            NO_n = u(idx.NO_n, :);
            
       end
        
        % Input of K+ into SC
        function f = input_K(self, t)      
            p = self.params;                
            f = zeros(size(t));
            if p.KSwitch == 1
                ii = p.t_0 <= t & t < p.t_1;
                f(ii) = ...
                    p.F_input * p.gab / ...
                    (p.ga * p.gb) * ...
                    (1 - (t(ii) - p.t_0) / p.delta_t).^(p.beta - 1) .* ...
                    ((t(ii) - p.t_0) / p.delta_t).^(p.alpha - 1);
                f(p.t_2 <= t & t <= p.t_3) = -p.F_input;
            end
        end
        
        function Glu = input_Glu(self, t) 
            p = self.params;
            Glu = p.GluSwitch * (p.Glu_max - p.Glu_min) * ( ...
            0.5 * tanh((t - p.t_0_Glu) / p.theta_L_Glu) - ...
            0.5 * tanh((t - p.t_2_Glu) / p.theta_R_Glu)) + p.Glu_min;
        end
        
        function names = varnames(self)
            names = [fieldnames(self.index); fieldnames(self.idx_out)];
        end
    end
end    
        
function idx = indices()    %for state variables
    idx.Ca_n = 1;                  
    idx.nNOS_act_n = 2;
    idx.NO_n = 3;
end        

function [idx, n] = output_indices()    %for variables in nargout loop
    idx.ft = 1;
    idx.J_Na_n = 2;
    idx.Glu = 3; 
    idx.w_NR2A = 4; 
    idx.w_NR2B = 5; 
    idx.I_Ca = 6; 
    idx.I_Ca_tot = 7; 
    idx.CaM = 10; 
    idx.tau_nk = 11; 
    idx.p_NO_n = 12;    
    idx.c_NO_n = 13; 
    idx.d_NO_n = 14;       
    
    n = numel(fieldnames(idx));
end
        
function params = parse_inputs(varargin)
    parser = inputParser();
    
    parser.addParameter('GluSwitch', 1); 
    parser.addParameter('KSwitch', 1); 
    parser.addParameter('NOswitch', 1); 
    
    % global constants
    parser.addParameter('F', 9.65e4); %C mol^-1; Faraday's constant
    parser.addParameter('R_gas', 8.315); %J mol^-1 K^-1; Gas constant
    parser.addParameter('T', 300); % K; Temperature
    
    % input 
    parser.addParameter('startpulse', 200);    
    parser.addParameter('lengthpulse', 200);
    parser.addParameter('lengtht1', 10);
    parser.addParameter('F_input', 2.5);        %s  
    parser.addParameter('alpha', 2);
    parser.addParameter('beta', 5);
    parser.addParameter('delta_t', 10);         %s
    parser.addParameter('k_C', 7.35e-5);        %uM m s^-1
    parser.addParameter('Glu_max', 1846);       % microM (one vesicle, Santucci2008)
    parser.addParameter('Glu_min', 0);          % microM
    parser.addParameter('theta_L_Glu', 1);      % slope of Glu input 
    parser.addParameter('theta_R_Glu', 1);      % slope of Glu input 

    % NO pathway 
    parser.addParameter('m_c', 4);              % [-] Number of Ca2+ bound per calmodulin (approximated as parameter, originally an algebraic variable)
    
    parser.addParameter('K_mA', 650);           % [uM] - fit to Santucci2008
    parser.addParameter('K_mB', 2800);          % [uM] - fit to Santucci2008

    parser.addParameter('v_n', -0.04);          % [V] ; the neuronal membrane potential , assumed to be approx constant in this model
    parser.addParameter('G_M', 46000);          % [fS]! was 46 pS! ; the conductance of the NMDA channel to Ca2+ compaired  
    parser.addParameter('P_Ca_P_M', 3.6);       % [-] ; the relative conductance of the NMDA channel to Ca2+ compared to monovalent ions
    parser.addParameter('Ca_ex', 2e3);          % [microM] ; the external calcium concentration (in Comerford+David2008: 1.5 mM!)
    parser.addParameter('M', 1.3e5);            % [microM] ; the concentration of monovalent ions in the neuron

    parser.addParameter('n_NR2A', 0.63);        % [-] ; average number of NR2A NMDA receptors per synapse (Santucci2008)
    parser.addParameter('n_NR2B', 11);          % [-] ; average number of NR2B NMDA receptors per synapse (Santucci2008)
    parser.addParameter('V_max_NO_n', 4.22);    % [s^-1] ; maximum catalytic rate of NO production (Chen2006) - obtained from fig 6 & equ 17 & 18
    parser.addParameter('O2_n', 200);           % [uM] ; tissue O2 concentration in the neuron (M.E.)
    parser.addParameter('K_mO2_n', 243);        % [uM] ; Chen2006
    parser.addParameter('LArg_n', 100);         % [uM] ; 
    parser.addParameter('K_mArg_n', 1.5);       % [uM] ; 
    
    parser.addParameter('k_O2_n', 9.6e-6);      % [uM^-2 s^-1] ; % (Kavdia2002)
    parser.addParameter('x_nk', 25);            % [um] ;  (M.E.)
    parser.addParameter('D_cNO', 3300);         % [um^2 s^-1] ; Diffusion coefficient NO (Malinski1993)
    
    parser.addParameter('V_spine', 8e-8);       % [nL] ; volume of the neuronal dendritic spine Santucci2008
    parser.addParameter('k_ex', 1600);          % [s^-1] ; decay rate constant of internal calcium concentration Santucci2008
    parser.addParameter('Ca_rest', 0.1);        % [uM] ; resting calcium concentration (in Comerford+David2008: 2.830 mM; in Santucci2008P: 0.1 \muM)
    parser.addParameter('lambda_buf', 20);      % [-] ; buffer capacity Santucci2008

    parser.addParameter('V_maxNOS', 25e-3);     % [] ; M.E.
    parser.addParameter('K_actNOS', 9.27e-2);   % [uM] ; 
    parser.addParameter('mu2_n', 0.0167);       % [s^-1] ; rate constant at which the nNOS is deactivated Comerford2008
    
    parser.addParameter('Q1', 1.9e5);           % [uM^-1]
    parser.addParameter('Q2', 2.1e5);           % [uM^-1]
    parser.addParameter('Q3', 0.4e5);           % [uM^-1]
    parser.addParameter('Q4', 0.26e5);          % [uM^-1]



    parser.parse(varargin{:})
    params = parser.Results;
    params.t_0 = params.startpulse;
    params.t_1 = params.t_0 + params.lengtht1;
    params.t_2 = params.t_0 + params.lengthpulse;
    params.t_3 = params.t_1 + params.lengthpulse;
    params.gab = factorial(params.alpha + params.beta - 1);
    params.ga = factorial(params.alpha - 1);
    params.gb = factorial(params.beta - 1);
    params.t_0_Glu = params.t_0;
    params.t_2_Glu = params.t_2;
    params.startpulse_Glu = params.startpulse;
    params.lengthpulse_Glu = params.lengthpulse;
end

function u0 = initial_conditions(idx)
    u0 = zeros(length(fieldnames(idx)), 1);
    
    u0(idx.Ca_n) = 0.0001;
    u0(idx.nNOS_act_n) = 0.3;
    u0(idx.NO_n) = 0.1;
end
