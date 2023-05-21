clear
close all
%Choose model P0, P1, 2RC, 3RC, 3RC-P0
model = "P1";

%Choose the sample cell(s)
cells = 1:3;

%Choose the starting point(s)
start_points = 1:2;

%Choose the initial offset on the SoC.
% -1 -> -5%
% 0 -> 0%, correct
% 1 -> +5%
errors_soc = -1:1;

% Do not change this parameter, use the file EKF_Marokko
marokko = false;

%Choose if visuals should be displayed
visuals = true;

results = zeros(length(cells), length(errors_soc), length(start_points));

for start_point = start_points
    for error_soc = errors_soc
        if start_point == 1
            offset_index = 11956;
            est_current = 0;
            est_current3 = 500;
            est_soc = 0.3996;
        elseif start_point == 2
            est_soc = 0.5448;
            offset_index = 19511;
            est_current = 1000;
            est_current3 = 1000;
        end
    
        est_soc = est_soc + error_soc * 0.05;

        Rk = 1e5;
        rho = 1e1;
        
        std_soc = 0.05;
        std_current = 1000;
        
        est_volt = 0;
        std_volt = 1;

        for cell = cells
            model
            start_point
            error_soc
            cell
            str_input = strcat("../../../Data/RACE/Cel", num2str(cell), "/input_race_full.mat");
            input = struct2cell(load(str_input));
            input = input{1};
            UKF_main;
        end
    end
end
