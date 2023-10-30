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

% ---------- DELTA SIRD Model ----------

% Needed for lsim
B = zeros(4,1); 

% ----- Initial guesses -----
% Rates
ki = 0.05; % Infection rate
kr = 0.5; % Recovery rate
kd = 0.01; % Death rate

% Initial conditions
S0 = .9;
I0 = .05;
R0 = .05;
D0 = 1-(S0+I0+R0);

deltaStart = [ki kr kd S0 I0 R0];

% ----- Constraints -----
lb = zeros(1,length(deltaStart)); % Lower bound
ub = ones(1,length(deltaStart)); % Upper bound


Aeq = [%0 0 0 -1 -1 -1; % -S0 - I0 - R0 < 2
       %0 0 0 1 1 1     %  S0 + I0 + R0 < 1
       0 1 1 0 0 0];    %  kr + kd < 1
beq = [%2
       %1
        1];

% Aeq = [0 1 1 0 0 0];
% beq = [1];

% ----- Optimization -----

% Capture local variables and allow multiple lines in objective function
deltaClosure = @(deltaA) objective_function(deltaA, delta_cases_prop, delta_deaths_prop);

% Display fmincon progress
% options = optimoptions('fmincon', 'Display', 'iter');

% Perform the optimization
delta_optimal = fmincon(deltaClosure, deltaStart, [], [], Aeq, beq, lb, []);

% ----- Results -----

% Extract optimization results
ki = delta_optimal(1);
kr = delta_optimal(2);
kd = delta_optimal(3);
S0 = delta_optimal(4);
I0 = delta_optimal(5);
R0 = delta_optimal(6);

% Create LTI
A = [1-ki 0 0 0;
     ki 1-(kd+kr) 0 0;
     0 kr 1 0;
     0 kd 0 1];
initialState = [S0 I0 R0 1-(S0+I0+R0)];

% Run the optimized model
delta_sys_sir_base = ss(A,B,eye(4),zeros(4,1),1);
deltaTotalModel = lsim(delta_sys_sir_base,zeros(16,1),linspace(0,15,16),initialState);

% Convert to form of real world data
deltaModeldCumulativeCases = (1-S0)+cumsum(ki*deltaTotalModel(:,1))';

figure
plot(delta_dates, deltaModeldCumulativeCases);
hold on
plot(delta_dates, delta_cases_prop);
title("Modeled vs. Actual Covid Cases")
legend("Model", "Actual")
xlabel("Date")
ylabel("Total Cases")

figure
plot(delta_dates, deltaTotalModel(:,4));
hold on
plot(delta_dates, delta_deaths_prop);
legend("Model", "Actual")



function error = objective_function(input, actualCases, actualDeaths) 
    % Extract and name inputs
    ki = input(1);
    kr = input(2);
    kd = input(3);
    S0 = input(4);
    I0 = input(5);
    R0 = input(6);
    D0 = 1-(S0+I0+R0);
    
    % Create LTI
    A = [1-ki 0 0 0;
         ki 1-(kd+kr) 0 0;
         0 kr 1 0;
         0 kd 0 1];
    initialState = [S0 I0 R0 D0];

    % Run SIRD model with current parameters
    sys_sir_base = ss(A, zeros(4,1), eye(4), zeros(4,1),1);
    modeled = lsim(sys_sir_base, zeros(16,1),linspace(0,16-1,16),initialState);

    % Convert to cumulative cases
    modeledCumulativeCases = (1-S0)+cumsum(ki*(modeled(:,1)));

    % Compare to real world data
    infectionError = norm(modeledCumulativeCases - actualCases);
    deathError = norm(modeled(:,4) - actualDeaths);

    % Report error
    error = infectionError+deathError;
end