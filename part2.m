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