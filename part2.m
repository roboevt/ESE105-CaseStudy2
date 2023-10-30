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

delta_deaths = deaths_STL(delta_idx);
omicron_deaths = deaths_STL(omicron_idx);

% Convert to proportions
delta_cases_prop = delta_cases/POP_STL;
delta_deaths_prop = delta_deaths/POP_STL;
omicron_cases_prop = omicron_cases/POP_STL;
omicron_deaths_prol = omicron_deaths/POP_STL;

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

% DELTA SIRD Model

deltaA = [0.95 0.04 0.01 0; 
    0.05 0.85 0.01 0; 
    0 0.1 0.98 0; 
    0 0.01 0 1];

delta0 = [1-delta_prop(1) delta_prop(1) 0 0];

% Needed for lsim
deltaB = zeros(4,1); 

% Capture local variables and allow multiple lines in objective function
deltaClosure = @(deltaA) objective_function(deltaA, deltaB, delta_cases_prop, delta_deaths_prop, delta0);

% Define the lower and upper bounds for A
lb = ones(4, 4)*0.00000;
ub = ones(4, 4)*0.99999;

% Ensure that the sum of every row and column of A is 1

AeqHelper = [1 0 0 0 1 0 0 0 1 0 0 0 1]; %Extract a row

Aeq = [ones(1,4) zeros(1,12);           %First column
       zeros(1,4) ones(1,4) zeros(1,8); %Second column
       zeros(1,8) ones(1,4) zeros(1,4); %Third Column
       zeros(1,12) ones(1,4);           % Fourth Column
       AeqHelper 0 0 0;    %First row
       0 AeqHelper 0 0;    %Second row
       0 0 AeqHelper 0;    %Third row
       0 0 0 AeqHelper; ]; %Fourth row

Beq = ones(1,16); %All sum to 1

% Display fmincon progress
% options = optimoptions('fmincon', 'Display', 'iter');

% Perform the optimization
delta_optimal = fmincon(deltaClosure, deltaA, [], [], Aeq, beq, lb, ub, [], options);

delta_sys_sir_base = ss(delta_optimal,deltaB,eye(4),zeros(4,1),1);
deltaTotalModel = lsim(delta_sys_sir_base,zeros(16,1),linspace(0,15,16),delta0);
deltaModeldCumulative = cumsum(deltaTotalModel(:,2)*POP_STL)';

figure
plot(delta_dates, deltaModeldCumulative);
hold on
plot(delta_dates, delta_cases)
title("Modeled vs. Actual Covid Cases")
legend("Model", "Actual")
xlabel("Date")
ylabel("Total Cases")



function error = objective_function(A, B, actualCases, actualDeaths, delta0) 
    sys_sir_base = ss(A, B, eye(4), zeros(4,1),1);
    modeled = lsim(sys_sir_base, zeros(16,1),linspace(0,16-1,16),delta0);
    modeledCumulativeCases = cumsum(modeled(:,2))';
    modeledCumulativeDeaths = cumsum(modeled(:,4))';
    infectionSum = sum((modeledCumulativeCases-actualCases).^2);
    deathSum = sum((modeledCumulativeDeaths-actualDeaths).^2);
    error = infectionSum+deathSum;
end