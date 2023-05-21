%% P1
clear
close all

parameters = struct2cell(load("Parameters/P1/Best_P1-19_03-09_04.mat"));

% Choose the input-output sequence to compare the model with
input = struct2cell(load("../Data/RACE/Cel1/input_race_full.mat"));
%input = struct2cell(load("../Data/RACE/Cel2/input_race_full.mat"));
%input = struct2cell(load("../Data/RACE/Cel3/input_race_full.mat"));

parameters = parameters{1};
SOC_levels = parameters.("SOC_levels");
%SOC_OCV = SOC_levels;
OCV = parameters.("OCV");
gamma_ = parameters.("gamma");
R0 = parameters.("R0");
M0 = zeros(height(parameters), 1);
M1 = parameters.("M1");
M2 = parameters.("M2");
M3 = parameters.("M3");
M4 = parameters.("M4");
hyst_M = [M0 M1 M2 M3 M4];
current_hyst = [0 750 3750 7500 15000];

model = "Models/Model_P1.slx";

%%
input = input{1};
input.("Voltage(V)") = input.("Voltage(V)")*1000;

I_rc = input.("Current(mA)");
I_rc = timeseries(I_rc, input.Time);
SOC_ts = input.("SOC");
SOC_ts = timeseries(SOC_ts, input.Time);
duration_vector = timeseries([0; diff(input.Time)], input.Time);

out_final = sim(model, max(I_rc.Time));
%out = sim(model, 20);

V0_final = out_final.yout{1}.Values;
hyst_final = out_final.yout{4}.Values;
result_final = out_final.yout{2}.Values;
ocv_final = out_final.yout{3}.Values;

result_final = timeseries(interp1(result_final.Time, result_final.Data, input.Time, "linear", "extrap"), input.Time);
ocv_final = timeseries(interp1(ocv_final.Time, ocv_final.Data, input.Time, "linear", "extrap"), input.Time);
V0_final = timeseries(interp1(V0_final.Time, V0_final.Data, input.Time, "linear", "extrap"), input.Time);
hyst_final = timeseries(interp1(hyst_final.Time, hyst_final.Data, input.Time, "linear", "extrap"), input.Time);

start_race = find(I_rc.Data(11000:end) > 0, 1,'first') + 11000;
start_hppc = find(I_rc.Data(23001:end) < 0, 1,'first') + 23000;

error = timeseries(result_final.Data - input.("Voltage(V)"), input.Time);
error_abs = timeseries(abs(error.Data), error.Time);
mov_error = timeseries(movmean(error.Data, 60), error.Time);
error_2 = timeseries(error_abs.Data.^2, error.Time);

error_hppc = timeseries(error_abs.Data(start_hppc:end, :), error_abs.Time(start_hppc:end, :));
error_2_hppc = timeseries(error_2.Data(start_hppc:end, :), error_2.Time(start_hppc:end, :));
error_race = timeseries(error_abs.Data(start_race:end, :), error_abs.Time(start_race:end, :));
error_2_race = timeseries(error_2.Data(start_race:end, :), error_2.Time(start_race:end, :));
result_hppc = timeseries(result_final.Data(start_hppc:end, :), result_final.Time(start_hppc:end, :));
result_race = timeseries(result_final.Data(start_race:end, :), result_final.Time(start_race:end, :));

duration_vector = [0; diff(I_rc.Time)];
duration_vector_hppc = duration_vector(start_hppc:end, :);
duration_vector_race = duration_vector(start_race:end, :);

RMSE = sqrt(sum(error_2.Data .* duration_vector)/max(error_2.Time));
ME = sum(error_abs.Data .* duration_vector)/max(error_abs.Time);
RMSE_RACE = sqrt(sum(error_2_race.Data .* duration_vector_race)/max(error_2.Time));
RMSE_HPPC = sqrt(sum(error_2_hppc.Data .* duration_vector_hppc)/max(error_2.Time));
ME_RACE = sum(error_race.Data .* duration_vector_race)/max(error_abs.Time);
ME_HPPC = sum(error_hppc.Data .* duration_vector_hppc)/max(error_abs.Time);

display("ME: " + ME + " mV");
display("RMSE: " + +RMSE + " mV");
display("ME without preprocessing: " + ME_RACE + " mV");
display("RMSE without preprocessing: " + RMSE_RACE + " mV");
display("Error end of Day: " + error.Data(end) + " mV");

figure()
sgtitle("P1-model")
ax1 = subplot(4, 1, 1);
plot(result_final.Time, result_final.Data);
hold on
plot(input.Time, input.("Voltage(V)"))
hold on
plot(ocv_final.Time, ocv_final.Data)
grid()
ylabel("Voltage [mV]")
title("Model vs. Target Voltage")
xlabel("Time [s]")
legend("Model", "Target", "OCV")
ax2 = subplot(4,1,2);
plot(V0_final.Time, V0_final.Data)
hold on
plot(hyst_final.Time, hyst_final.Data)
legend("R0", "Delta")
ylabel("Voltage [mV]")
title("Voltage per impedance")
xlabel("Time [s]")
grid()
ax3 = subplot(4, 1, 3);
plot(error.Time, error.Data);
hold on
plot(mov_error.Time, mov_error.Data)
grid()
legend("Error", "Moving Error")
ylabel("Voltage [mV]")
title("Error Voltage")
xlabel("Time [s]")
ax4 = subplot(4, 1, 4);
plot(SOC_ts.Time, SOC_ts.Data)
linkaxes([ax1, ax2, ax3, ax4], 'x')
grid()
title("SoC")
ylabel("SoC [-]")
xlabel("Time [s]")
grid("on")

figure()
plot(result_final.Time, result_final.Data, "LineWidth", 1);
hold on
plot(input.Time, input.("Voltage(V)"), "LineWidth", 1)
hold on
plot(result_final.Time, ocv_final.Data, "LineWidth", 1)
xlabel("Time [s]")
ylabel("Voltage [mV]")
title("Model vs. Target Voltage")
legend("Model", "Target", "OCV")
grid("on")

