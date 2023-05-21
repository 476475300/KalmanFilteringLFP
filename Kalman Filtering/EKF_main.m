%% Framework for Extended Kalman Filter
%%
dir = "../Parameters/";
folder = strcat("../Parameters/", model);
if model == "2RC"
    parameters = struct2cell(load(strcat(folder, "/Best_2RC_03_03-15_44.mat")));    
elseif model == "3RC"
    parameters = struct2cell(load(strcat(folder, "/Best_3RC-03_03-22_05.mat")));
elseif model == "P0"
    parameters = struct2cell(load(strcat(folder, "/Best_P0_11_03-08_56.mat")));
elseif model == "P1"
    parameters = struct2cell(load(strcat(folder, "/Best_P1-19_03-09_04.mat")));
elseif model == "3RC-P0"
    parameters = struct2cell(load(strcat(folder, "/Best_3RC_P0-01_04-21_31")));
end
ocv_table = strcat(dir, "OCV_SOC");
ocv_table = table_from_csv(ocv_table);

parameters = parameters{1};
parameters.OCV = ocv_table.Avg_avg_OCV*1000;
parameters.dOCV = ocv_table.Avg_avg_dOCV*1000*1000;

if marokko
    input.("Voltage(V)") = input.(strcat("V_Cell", num2str(cell)))*1000;
else
    input.("Voltage(V)") = input.("Voltage(V)")*1000;
end

input = input(offset_index:end, :);
time_steps = [0; diff(input.Time)];

%% Initialisation
SOC_real = [input.("SOC")];
est_soc = SOC_real(1) + error_soc*0.05;
init_error = est_soc - SOC_real(1);
SOC_Uncertainty = zeros(height(input), 1);

if model == "2RC"
    Pk_k = diag([std_soc^2, std_current^2, std_current^2]); 
    xk_k = [est_soc, est_current, est_current]';
elseif model == "3RC"
    Pk_k = diag([std_soc^2, std_current^2, std_current^2, std_current^2]); 
    xk_k = [est_soc, est_current, est_current, est_current3]';
elseif model == "P0"
    Pk_k = diag([std_soc^2]); 
    xk_k = [est_soc]';
elseif model == "P1"
    Pk_k = diag([std_soc^2 std_volt^2]); 
    xk_k = [est_soc est_volt]';
elseif model == "Comb0"
    Pk_k = diag([std_soc^2, std_current^2, std_current^2]); 
    xk_k = [est_soc, est_current, est_current]';
elseif model == "3RC-P0"
    Pk_k = diag([std_soc^2, std_current^2, std_current^2, std_current^2]); 
    xk_k = [est_soc, est_current, est_current, est_current3]';
end
SOC_kalman = [xk_k(1); zeros(height(input)-1, 1)];
n = length(xk_k);
%%
V_kalman = zeros(height(input), 1);
V_real = [input.("Voltage(V)")];

if model == "P0" || model == "Comb0" || model == "3RC-P0"
    V_rc = zeros(height(input), n+2);
    state_hyst = 1;
else
    V_rc = zeros(height(input), n+1);
end

%2D interpolation of SOC and current for M
if model == "P1"
    current_hyst = [0; 1500; 3750; 7500; 15000];
    sampling_current_hyst = 0:100:15000;
    M_matrix = interp2(current_hyst, parameters.("SOC_levels"), [zeros(height(parameters), 1) parameters{:, 6:end}], sampling_current_hyst, parameters.("SOC_levels"),"linear");
end

OCV_exp = parameters.("OCV")(round(input.("SOC")*1000)+1);

w = waitbar(0, "Starting");
%output_function = eval(output_function);
for i = 1:height(input)
    if rem(i, 1000) == 0
        waitbar(i/height(input), w, "Processing Cell " + num2str(cell) +  " (" + round(i/height(input), 2)*100 + "%)")
    end
    current = input.("Current(mA)")(i);

    if model == "2RC"
        [A, B, C, R0, R1, R2, OCV] = load_matrices(cell, model, parameters, xk_k(1), current, time_steps(i), marokko);
    elseif model == "3RC"
        [A, B, C, R0, R1, R2, R3, OCV] = load_matrices(cell, model, parameters, xk_k(1), current, time_steps(i), marokko);
    elseif model == "P0"
        [A, B, C, OCV_Delta, R0, epsilon, OCV] = load_matrices(cell, model, parameters, xk_k(1), current, time_steps(i), marokko);
        if current > epsilon
            state_hyst = 1;
        elseif current < -1*epsilon
            state_hyst = -1;
        end
    elseif model == "P1"
        [A, B, C, OCV, R0, gamma, M] = load_matrices(cell, model, parameters, xk_k(1), current,  time_steps(i), marokko, xk_k(2), M_matrix, sampling_current_hyst);
    elseif model == "Comb0"
        [A, B, C, R0, R1, R2, OCV, OCV_offset, epsilon] = load_matrices(cell, model, parameters, xk_k(1), current, time_steps(i), marokko);
        if current > epsilon
            state_hyst  = 1;
        elseif current < -1*epsilon
            state_hyst = -1;
        end
    elseif model == "3RC-P0"
        [A, B, C, R0, R1, R2, R3, OCV, OCV_offset, epsilon] = load_matrices(cell, model, parameters, xk_k(1), current, time_steps(i), marokko);
        if current > epsilon
            state_hyst  = 1;
        elseif current < -1*epsilon
            state_hyst = -1;
        end
    end
    Qk = B*rho^2*B';

    % Predictor
    if model == "P1"
        xk_kmin = state_function(xk_k, [current; M], A, B);
    else
        xk_kmin = state_function(xk_k, current, A, B);
    end
    Pk_kmin = A*Pk_k*A' + Qk;

    % Predict output model
    if model == "2RC"
        yk = output_function(model, current, OCV, xk_kmin, R0, R1, R2);
    elseif model == "3RC"
        yk = output_function(model, current, OCV, xk_kmin, R0, R1, R2, R3);
    elseif model == "P0"
        yk = output_function(model, current, OCV, OCV_Delta, R0, state_hyst);   
    elseif model == "P1"
        yk = output_function(model, current, OCV, xk_kmin, R0);
    elseif model == "Comb0"
        yk = output_function(model, current, OCV, xk_kmin, R0, R1, R2, OCV_offset, state_hyst);
    elseif model == "3RC-P0"
        yk = output_function(model, current, OCV, xk_kmin, R0, R1, R2, R3, OCV_offset, state_hyst);       
    end

    % Calculate error
    yk_real = V_real(i);
    yk_error = yk_real - yk;

    % Update
    Sk = C*Pk_kmin*C' + Rk;
    Kk = Pk_kmin*C'*Sk^(-1);
    xk_k = xk_kmin + Kk*yk_error;
    Pk_k = (eye(n) - Kk*C)*Pk_kmin;

    % Collection vectors
    V_kalman(i) = yk;
    SOC_kalman(i) = xk_k(1);

    if model == "2RC"
        voltages = xk_kmin(2:end) .* [R1; R2];
        V_rc(i, :) = [OCV current*R0 voltages'];
    elseif model == "3RC"
        voltages = xk_kmin(2:end) .* [R1; R2; R3];
        V_rc(i, :) = [OCV current*R0 voltages'];
    elseif model == "P0"
        V_rc(i, :) = [OCV current*R0 state_hyst*OCV_Delta];
    elseif model == "P1"
        V_rc(i, :) = [OCV current*R0 xk_k(2)];
    elseif model == "Comb0"
        voltages = xk_kmin(2:end) .* [R1; R2];
        V_rc(i, :) = [OCV current*R0 voltages' state_hyst*OCV_offset];
    elseif model == "3RC-P0"
        voltages = xk_kmin(2:end) .* [R1; R2; R3];
        V_rc(i, :) = [OCV current*R0 voltages' state_hyst*OCV_offset];
    end
    SOC_Uncertainty(i) = sqrt(Pk_k(1, 1));
end
close(w)
%%
SOC_error = SOC_kalman - SOC_real;
V_error = V_kalman - V_real;
output = table(SOC_kalman, SOC_real, SOC_error, SOC_Uncertainty, V_kalman, V_real, V_error);
%%
if visuals
    figure
    super_title = strcat("Performance UKF: ", model);
    sgtitle(super_title);
    ax1 = subplot(4, 1, 1);
    plot(input.Time, output.("V_real"), "LineWidth", 1)
    hold on
    plot(input.Time, output.("V_kalman"), "LineWidth", 1)
    legend("Real", "Kalman")
    xlabel("Time [s]")
    ylabel("Voltage [mV]")
    title("Kalman vs. Target voltage")
    grid on
    ax2 = subplot(4, 1, 2);
    plot(input.Time, output.("SOC_real"), "LineWidth", 1)
    hold on
    plot(input.Time, output.("SOC_kalman"), "LineWidth", 1)
    legend("Real", "Kalman")
    title("SoC")
    xlabel("Time [s]")
    ylabel("SoC [-]")
    title("Kalman vs. Target SoC")
    grid on
    ax3 = subplot(4, 1, 3);
    plot(input.Time, output.("V_error"), "LineWidth", 1)
    title("Voltage Error")
    xlabel("Time [s]")
    ylabel("Voltage [mV]")
    title("Error on voltage")
    grid on
    ax4 = subplot(4, 1, 4);
    plot(input.Time, output.("SOC_error"), "LineWidth", 1)
    title("SoC Error")
    linkaxes([ax1 ax2 ax3 ax4], 'x')
    xlabel("Time [s]")
    ylabel("SoC [-]")
    title("Error on SoC")    
    grid on
    
    figure()
    super_title = strcat("Output UKF: ", model);
    sgtitle(super_title)
    ax1 = subplot(2, 1, 1);
    plot(input.Time, output.("V_real"))
    hold on
    plot(input.Time, output.("V_kalman"))
    hold on
    plot(input.Time, V_rc(:, 1))
    hold on
    plot(input.Time, OCV_exp)
    legend("Real", "Kalman", "OCV", "Expected OCV")
    grid on
    title("Kalman vs. Target voltage")
    ax2 = subplot(2, 1, 2);
    for i = 2:size(V_rc, 2)
        plot(input.Time, V_rc(:, i))
        hold on
    end
    if model == "2RC"
        legend("R0", "RC1", "RC2")
    elseif model == "3RC"
        legend("R0", "RC1", "RC2", "RC3")
    elseif model == "P0"
        legend("R0", "Hyst")
    elseif model == "P1"
        legend("R0", "Hyst")
    elseif model == "Comb0"
        legend("R0", "RC1", "RC2", "Hyst")    
    elseif model == "3RC-P0"
        legend("R0", "RC1", "RC2", "RC3", "Hyst")    
    end    
    grid on
    title("Voltage per impedance")
    linkaxes([ax1 ax2], 'x')    
end
%%
error_avg = sum(abs(SOC_error) .* time_steps)/sum(time_steps);
error_2_avg = sqrt(sum((SOC_error.^2 .* time_steps)/sum(time_steps)));
display("Initial Error: " + init_error*100 + " %")
display("Avg Error: "+ error_avg * 100 + " %")
display("RMSE: " + error_2_avg * 100 + " %")
display("Final Error: " + SOC_error(end) * 100 + " %")
