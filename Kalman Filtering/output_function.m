function terminal_volt = output_function(model, current, OCV, varargin)
    if model == "2RC"
        state = varargin{1};
        R0 = varargin{2}; 
        R1 = varargin{3};
        R2 = varargin{4};
        terminal_volt = OCV + current*R0 + [R1 R2] * state(2:end);   
    elseif model == "3RC"
        state = varargin{1};
        R0 = varargin{2}; 
        R1 = varargin{3};
        R2 = varargin{4};
        R3 = varargin{5};
        terminal_volt = OCV + current*R0 + [R1 R2 R3] * state(2:end);
    elseif model == "P0"
        OCV_Delta = varargin{1};
        R0 = varargin{2};
        state_hyst = varargin{3};
        terminal_volt = OCV + current*R0 + state_hyst*OCV_Delta;
    elseif model == "P1"
        V_hyst = varargin{1}(2);
        R0 = varargin{2};
        terminal_volt = OCV + current*R0 + V_hyst;
    elseif model == "Comb0"
        state = varargin{1};
        R0 = varargin{2}; 
        R1 = varargin{3};
        R2 = varargin{4};
        OCV_offset = varargin{5};
        state_hyst = varargin{6};
        terminal_volt = OCV + current*R0 + [R1 R2] * state(2:end) + state_hyst*OCV_offset;   
    elseif model == "3RC-P0"
        state = varargin{1};
        R0 = varargin{2}; 
        R1 = varargin{3};
        R2 = varargin{4};
        R3 = varargin{5};
        OCV_offset = varargin{6};
        state_hyst = varargin{7};
        terminal_volt = OCV + current*R0 + [R1 R2 R3] * state(2:end) + state_hyst*OCV_offset;   
    end
end