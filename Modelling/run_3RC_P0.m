%% 3RC-P0
clear
close all

parameters = struct2cell(load("Parameters/3RC-P0/Best_3RC_P0-01_04-21_31.mat"));

% Choose the input-output sequence to compare the model with
input = struct2cell(load("../Data/RACE/Cel1/input_race_full.mat"));
%input = struct2cell(load("../Data/RACE/Cel2/input_race_full.mat"));
%input = struct2cell(load("../Data/RACE/Cel3/input_race_full.mat"));

% Parameters are loaded
parameters = parameters{1};
SOC_levels = parameters.("SOC_levels");
OCV = parameters.("OCV");
tau1 = parameters.("tau1");
tau2 = parameters.("tau2");
tau3 = parameters.("tau3");
R0 = parameters.("R0");
R1 = parameters.("R1");
R2 = parameters.("R2");
R3 = parameters.("R3");

epsilon = parameters.("epsilon");
OCV_offset = parameters.("OCV_offset");

model = "Models/Model_3RC_P0.slx";

%%
input = input{1};
input.("Voltage(V)") = input.("Voltage(V)")*1000;

I_rc = input.("Current(mA)");
I_rc = timeseries(I_rc, input.Time);
SOC_ts = input.("SOC");
SOC_ts = timeseries(SOC_ts, input.Time);

out = sim(model, max(I_rc.Time));
rc1_out = out.yout{1}.Values;
rc2_out = out.yout{2}.Values;
rc3_out = out.yout{3}.Values;
result = out.yout{7}.Values;
ocv_out = out.yout{4}.Values;
r0_out = out.yout{5}.Values;
ocv_offset_out = out.yout{6}.Values;

result_int = interp1(result.Time, result.Data, input.Time, "linear", "extrap");
ocv_int = interp1(ocv_out.Time, ocv_out.Data, input.Time, "linear", "extrap");

start_race = find(I_rc.Data(11000:end) > 0, 1,'first') + 11000;
start_hppc = find(I_rc.Data(23001:end) < 0, 1,'first') + 23000;

error = timeseries(result_int - input.("Voltage(V)"), input.Time);
error_abs = timeseries(abs(error.Data), error.Time);
mov_error = timeseries(movmean(error.Data, 60), error.Time);
error_2 = timeseries(error_abs.Data.^2, error.Time);

error_hppc = timeseries(error_abs.Data(start_hppc:end, :), error_abs.Time(start_hppc:end, :));
error_2_hppc = timeseries(error_2.Data(start_hppc:end, :), error_2.Time(start_hppc:end, :));
error_race = timeseries(error_abs.Data(start_race:end, :), error_abs.Time(start_race:end, :));
error_2_race = timeseries(error_2.Data(start_race:end, :), error_2.Time(start_race:end, :));
result_hppc = timeseries(result.Data(start_hppc:end, :), result.Time(start_hppc:end, :));
result_race = timeseries(result.Data(start_race:end, :), result.Time(start_race:end, :));

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
sgtitle("3RC-P0-model")
ax1 = subplot(4, 1, 1);
plot(input.Time, result_int);
hold on
plot(input.Time, input.("Voltage(V)"))
hold on
plot(input.Time, ocv_int)
grid()
ylabel("Voltage [mV]")
title("Model vs. Target Voltage")
xlabel("Time [s]")
legend("Model", "Target", "OCV")
ax2 = subplot(4,1,2);
plot(r0_out.Time, r0_out.Data)
hold on
plot(rc1_out.Time, rc1_out.Data)
hold on
plot(rc2_out.Time, rc2_out.Data)
hold on
plot(rc3_out.Time, rc3_out.Data)
hold on
plot(ocv_offset_out.Time, ocv_offset_out.Data)
legend("RC1", "RC2", "RC3", "R0", "Offset")
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
plot(input.Time, result_int, "LineWidth", 1);
hold on
plot(input.Time, input.("Voltage(V)"), "LineWidth", 1)
hold on
plot(input.Time, ocv_int, "LineWidth", 1)
xlabel("Time [s]")
ylabel("Voltage [mV]")
title("Model vs. Target Voltage 3RC-P0")
legend("Model", "Target", "OCV")
grid("on")