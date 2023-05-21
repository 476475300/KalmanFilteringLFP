function [varargout] = lookup_values(model, parameters, SOC, varargin)
    [~, line_SOC] = min(abs(SOC-parameters.("SOC_levels")));
    OCV = parameters.("OCV")(line_SOC);
    dOCV = parameters.("dOCV")(line_SOC);
    if model == "2RC"
        varargout{1} = OCV;
        varargout{2} = dOCV;
        varargout{3} = parameters.("R0")(line_SOC);
        varargout{4} = parameters.("R1")(line_SOC);
        varargout{5} = parameters.("R2")(line_SOC);
        varargout{6} = parameters.("tau1")(line_SOC);
        varargout{7} = parameters.("tau2")(line_SOC);
    elseif model == "3RC"
        varargout{1} = OCV;
        varargout{2} = dOCV;
        varargout{3} = parameters.("R0")(line_SOC);
        varargout{4} = parameters.("R1")(line_SOC);
        varargout{5} = parameters.("R2")(line_SOC);
        varargout{6} = parameters.("R3")(line_SOC);
        varargout{7} = parameters.("tau1")(line_SOC);
        varargout{8} = parameters.("tau2")(line_SOC);
        varargout{9} = parameters.("tau3")(line_SOC);
    elseif model == "P0"
        varargout{1} = OCV;
        varargout{2} = parameters.("OCV_Delta")(line_SOC);
        varargout{3} = dOCV;
        varargout{4} = parameters.("R0")(line_SOC);
        varargout{5} = parameters.("epsilon")(line_SOC);
    elseif model == "P1"
        current = varargin{1};
        M_matrix = varargin{2};
        sampling_current_hyst = varargin{3};
        [~, line_current] = min(abs(current-sampling_current_hyst));
        M = M_matrix(line_SOC, line_current);
        varargout{1} = OCV;
        varargout{2} = parameters.("R0")(line_SOC);
        varargout{3} = dOCV;
        varargout{4} = parameters.("gamma")(line_SOC);
        varargout{5} = M;
    elseif model == "Comb0"
        varargout{1} = OCV;
        varargout{2} = dOCV;
        varargout{3} = parameters.("R0")(line_SOC);
        varargout{4} = parameters.("R1")(line_SOC);
        varargout{5} = parameters.("R2")(line_SOC);
        varargout{6} = parameters.("tau1")(line_SOC);
        varargout{7} = parameters.("tau2")(line_SOC);
        varargout{8} = parameters.("OCV_offset")(line_SOC);
        varargout{9} = parameters.("epsilon")(line_SOC);
    elseif model == "3RC-P0"
        varargout{1} = OCV;
        varargout{2} = dOCV;
        varargout{3} = parameters.("R0")(line_SOC);
        varargout{4} = parameters.("R1")(line_SOC);
        varargout{5} = parameters.("R2")(line_SOC);
        varargout{6} = parameters.("R3")(line_SOC);
        varargout{7} = parameters.("tau1")(line_SOC);
        varargout{8} = parameters.("tau2")(line_SOC);
        varargout{9} = parameters.("tau3")(line_SOC);
        varargout{10} = parameters.("OCV_offset")(line_SOC);
        varargout{11} = parameters.("epsilon")(line_SOC);
    end

end