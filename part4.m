close all

load("mockdata2023.mat")

% Font parameters
default_size = 18;
title_size = 24;
line_size = 3;

% Plot mockdata for reference, will add modeled data later
plot(newInfections, "LineWidth",line_size)
hold on
plot(cumulativeDeaths, "LineWidth", line_size)
title("Mock vs. Modeled Data", 'FontSize', title_size)
xlabel("Day", 'FontSize', default_size)
ylabel("Fraction of Population", 'FontSize', default_size)

% ---------- SIRVD Optimization ----------

% ----- Initial Guess -----

% Rates
ki = 0.002; % Infection rate
kr = 0.9939; % Recovery rate
kd = 0.061; % Death rate
kv = 0.01; % Vaccination rate
kb = 0.001; % Breakthrough infection rate
kp = 0.001; % ReInfection rate

% Start
S0 = 1; % Susceptible
I0 = 0; % Infected
R0 = 0; % Recovered
V0 = 0; % Vaccinated
D0 = 1 - S0 - I0 - R0 - V0; % Dead

d = 120; % Vaccine rollout day

initialGuess = [ki kr kd kv kb kp S0 I0 R0 V0 d];

% ----- Constraints -----
lb = zeros(1,length(initialGuess)); % Lower bound
ub = ones(1,length(initialGuess)); % Upper bound

Aeq = [0 0 0 0 0 0 1 1 1 1 0];
beq = [1];

% ----- Optimization -----

% Capture local variables and allow multiple lines in objective function
closure = @(state) objective_function(state, newInfections, cumulativeDeaths);

% Perform the optimization
optimal_result = fmincon(closure, initialGuess, [], [], Aeq, beq, lb, []);

% ---------- Optimization Results ----------
[modeledNewInfections, modeledCumulativeDeaths] = simulate(optimal_result);

plot(modeledNewInfections, "LineWidth",line_size)
plot(modeledCumulativeDeaths, "LineWidth", line_size)
legend("Mock New Infections", "Mock Cumulative Deaths", "Modeled New Infections", "Modeled Cumulative Deaths")

function [newInfections, cumulativeDeaths] = simulate(input)
    % Extract and name inputs
    ki = input(1);
    kr = input(2);
    kd = input(3);
    kv = 0; % We begin the simulation without vaccines
    kb = input(5);
    kp = input(6);
    S0 = input(7);
    I0 = input(8);
    R0 = input(9);
    V0 = input(10);
    D0 = 1-(S0+I0+R0+V0);
    d  = round(input(11));

    % Create LTI
    A = [1-ki-kv 0 kp 0 0;
        ki 1-(kd+kr) 0 kb 0;
        0 kr 1-(kp+kv) 0 0;
        kv 0 kv 1-kb 0;
        0 kd 0 0 1];
    initialState = [S0 I0 R0 V0 D0];
    
    % Run SIRD model with current parameters before vaccinations
    
    sys_sir_base = ss(A, zeros(5,1), eye(5), zeros(5,1),1);
    modeledBeforeVaccine = lsim(sys_sir_base, zeros(d,1),linspace(0,d-1,d),initialState);
    
    % Vaccine Rollout
    kv = input(4);
    A = [1-ki-kv 0 kp 0 0;
        ki 1-(kd+kr) 0 kb 0;
        0 kr 1-(kp+kv) 0 0;
        kv 0 kv 1-kb 0;
        0 kd 0 0 1];
    d = 400 - d;
    
    initialState(1) = modeledBeforeVaccine(1,end);
    initialState(2) = modeledBeforeVaccine(2,end);
    initialState(3) = modeledBeforeVaccine(3,end);
    initialState(4) = modeledBeforeVaccine(4,end);
    initialState(5) = modeledBeforeVaccine(5,end);
    
    sys_sir_base = ss(A, zeros(5,1), eye(5), zeros(5,1),1);
    modeledAfterVaccine = lsim(sys_sir_base, zeros(d,1),linspace(0,d-1,d),initialState);
    
    totalModel = [modeledBeforeVaccine; modeledAfterVaccine];
    
    % Convert results to match goal
    modeledNewInfections = totalModel(:,1) * ki;
    modeledBreakthroughInfections = totalModel(:,4) * kb;
    newInfections = modeledNewInfections + modeledBreakthroughInfections;
    
    cumulativeDeaths = totalModel(:,5);
end

% fmincon cost function
function error = objective_function(input, actualNewInfections, actualCumulativeDeaths)
    [totalModelNewInfections, modeledCumulativeDeaths] = simulate(input);
    
    % Compare to goal data
    infectionError = norm(totalModelNewInfections - actualNewInfections);
    deathError = norm(modeledCumulativeDeaths - actualCumulativeDeaths);
    
    % Report error
    error = infectionError+deathError;
end