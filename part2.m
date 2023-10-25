close all

load("COVID_STL.mat")

% Plot entire dataset of covid cases
plot(dates, cases_STL)
title("Cumulative Covid Cases in St. Louis")
xlabel("Date")
ylabel("Cases")

% Extract delta and omicron phases from data
delta_idx =  (dates > datetime("2021-6-30") & dates < datetime("2021-10-26"));
omicron_idx = (dates > datetime("2021-10-27") & dates < datetime("2022-3-22"));
delta_dates = dates(delta_idx);
omicron_dates = dates(omicron_idx);
delta_cases = cases_STL(delta_idx);
omicron_cases = cases_STL(omicron_idx);

% Plot delta and omicron phases for reference
figure
plot(delta_dates, delta_cases)
title("Delta Cases")
xlabel("Date")
ylabel("Cases")
figure
plot(omicron_dates, omicron_cases)
title("Omicron Cases")
xlabel("Date")
ylabel("Cases")

% SIRD Model

deltaA = [0.95 0.04 0.01 0; 
    0.05 0.85 0.01 0; 
    0 0.1 0.98 0; 
    0 0.01 0 1];

delta0 = [1 0 0 0];

sys_sir_base = ss(deltaA,B,eye(4),zeros(4,1),1);
Y = lsim(sys_sir_base,zeros(1000,1),linspace(0,999,1000),delta0);

figure
plot(Y)
legend('S','I','R','D');
xlabel('Time')
ylabel('Percentage Population');