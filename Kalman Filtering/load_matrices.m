function varargout = load_matrices(cell, model, parameters, current_soc, current, dt, marokko, varargin)
    if marokko
        C_nom = 15000;
    else
        C_nom = [14904, 15162, 15283];
        C_nom = C_nom(cell);
    end

    eta = 1;

    if model == "2RC"
        [OCV, dOCV, R0, R1, R2, tau1, tau2] = lookup_values(model, parameters, current_soc);
        A = diag([1 exp(-dt/tau1) exp(-dt/tau2)]);
        B = [eta*dt/(3600*C_nom); 1-exp(-dt/tau1); 1-exp(-dt/tau2)];
        C = [dOCV R1 R2];
        varargout{1} = A;
        varargout{2} = B;
        varargout{3} = C;
        varargout{4} = R0;
        varargout{5} = R1;
        varargout{6} = R2;
        varargout{7} = OCV;
    
    elseif model == "3RC"
        [OCV, dOCV, R0, R1, R2, R3, tau1, tau2, tau3] = lookup_values(model, parameters, current_soc);
        A = diag([1 exp(-dt/tau1) exp(-dt/tau2) exp(-dt/tau3)]);
        B = [eta*dt/(3600*C_nom); 1-exp(-dt/tau1); 1-exp(-dt/tau2); 1-exp(-dt/tau3)];
        C = [dOCV R1 R2 R3];
        varargout{1} = A;
        varargout{2} = B;
        varargout{3} = C;
        varargout{4} = R0;
        varargout{5} = R1;
        varargout{6} = R2;
        varargout{7} = R3;
        varargout{8} = OCV;

    elseif model == "P0"
        [OCV, OCV_Delta, dOCV, R0, epsilon]= lookup_values(model, parameters, current_soc);
        varargout{1} = 1;
        varargout{2} = eta*dt/(3600*C_nom);
        varargout{3} = dOCV;
        varargout{4} = OCV_Delta;
        varargout{5} = R0;
        varargout{6} = epsilon;
        varargout{7} = OCV;
        
    elseif model == "P1"
        current_hyst = varargin{1};
        M_matrix = varargin{2};
        sampling_current_hyst = varargin{3};
        [OCV, R0, dOCV, gamma, M]= lookup_values(model, parameters, current_soc, current, M_matrix, sampling_current_hyst);
        exp_term = exp(-1*abs(eta*current*gamma*dt/C_nom));
        varargout{1} = diag([1 exp_term]);
        varargout{2} = diag([eta*dt/(3600*C_nom) (1-exp_term)]);
        %varargout{3} = [dOCV gamma*sign(current)*(M-current_hyst)] ;
        varargout{3} = [dOCV 1];
        varargout{4} = OCV;
        varargout{5} = R0;        
        varargout{6} = gamma;
        varargout{7} = M;
    
    elseif model == "Comb0"
        [OCV, dOCV, R0, R1, R2, tau1, tau2, OCV_offset, epsilon] = lookup_values(model, parameters, current_soc);
        A = diag([1 exp(-dt/tau1) exp(-dt/tau2)]);
        B = [eta*dt/(3600*C_nom); 1-exp(-dt/tau1); 1-exp(-dt/tau2)];
        C = [dOCV R1 R2];
        varargout{1} = A;
        varargout{2} = B;
        varargout{3} = C;
        varargout{4} = R0;
        varargout{5} = R1;
        varargout{6} = R2;
        varargout{7} = OCV;
        varargout{8} = OCV_offset;
        varargout{9} = epsilon;
    elseif model == "3RC-P0"
        [OCV, dOCV, R0, R1, R2, R3, tau1, tau2, tau3, OCV_offset, epsilon] = lookup_values(model, parameters, current_soc);
        A = diag([1 exp(-dt/tau1) exp(-dt/tau2) exp(-dt/tau3)]);
        B = [eta*dt/(3600*C_nom); 1-exp(-dt/tau1); 1-exp(-dt/tau2); 1-exp(-dt/tau3)];
        C = [dOCV R1 R2 R3];
        varargout{1} = A;
        varargout{2} = B;
        varargout{3} = C;
        varargout{4} = R0;
        varargout{5} = R1;
        varargout{6} = R2;
        varargout{7} = R3;
        varargout{8} = OCV;
        varargout{9} = OCV_offset;
        varargout{10} = epsilon;
    end

end

