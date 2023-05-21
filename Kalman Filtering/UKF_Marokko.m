clear
close all

%Choose model P0, P1, 2RC, 3RC, 3RC-P0
model = "P1";

%Choose the module(s)
cells = 1:49;

%Choose the day
start_points = 1:5;

%Choose the initial offset on the SoC.
% -1 -> -5%
% 0 -> 0%, correct
% 1 -> +5%
error_soc = 0;

% Do not change this parameter, use the file EKF_script for lab data
marokko = true;

%Choose if visuals should be displayed
visuals = true;

soc_estimations = [0.993, 0.639, 0.477, 0.373, 0.301];  

for day = days  
    offset_index = 1;
    est_current = 0;
    est_current3 = 0;
    est_soc = soc_estimations(day);

    Rk = 1e5;
    rho = 1e1;
    
    std_soc = 0.05;
    std_current = 100;
    
    est_volt = 0;
    std_volt = 10;    
    str_input = strcat("../../../Data/MAROKKO/Clean_Data/Day", num2str(day), "_clean.mat");
    input = struct2cell(load(str_input));
    input = input{1};
    input.SOC = input.SOC_Ah/100;
    results_vector = zeros(height(input), 50);
    for cell = cells
        cell
        UKF_main;
        results_vector(:, cell) = SOC_kalman;
    end
    q1 = quantile(results_vector(:, cells), 0.25,2);
    q2 = quantile(results_vector(:, cells), 0.5,2);
    q3 = quantile(results_vector(:, cells), 0.75,2);
    q0 = min(results_vector(:, cells), [], 2);
    q4 = max(results_vector(:, cells),[], 2);

    results_vector(:, 50) = SOC_real;
    figure
    plot(input.Time, q0, "-.", "LineWidth", 1)
    hold on
    plot(input.Time, q1, "-.", "LineWidth", 1)
    plot(input.Time, q2, "LineWidth", 1.5)
    plot(input.Time, q3, "-.", "LineWidth", 1)
    plot(input.Time, q4, "-.", "LineWidth", 1)
    plot(input.Time, results_vector(:, 50), "LineWidth", 1.5)
    legend("Min", "Q1", "Median", "Q3", "Max", "Real")
    grid on
    xlabel("Time [s]")
    ylabel("SoC [-]")
    title(strcat("Accuracy of UKF & ", model, "-model  on data SCM: day ", num2str(day)))
end
%%



