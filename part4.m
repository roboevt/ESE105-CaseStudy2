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
title("Mock vs. Modeled SIRVD Data", 'FontSize', title_size)
xlabel("Day", 'FontSize', default_size)
ylabel("Fraction of Population", 'FontSize', default_size)

% ---------- SIRVD Optimization ----------

% ----- Initial Guess -----

% Rates
ki = 0.0058; % Infection rate
kr = 0.095; % Recovery rate
kd = .0136; % Death rate
kv = 0.04; % Vaccination rate
kb = 0.00002; % Breakthrough infection rate
kp = 0.001; % Reinfection rate

% Start point
S0 = 1; % Susceptible
I0 = 0; % Infected
R0 = 0; % Recovered
V0 = 0; % Vaccinated
D0 = 1 - S0 - I0 - R0 - V0; % Dead

d = 122; % Vaccine rollout day

initialGuess = [ki kr kd kv kb kp S0 I0 R0 d];

% ----- Constraints -----
lb = zeros(1,length(initialGuess)); % Lower bound
ub = ones(1,length(initialGuess)); % Upper bound
ub(end) = 400;
ub(9) = 0.1; % R0

Aeq = [0 0 0 0 0 0 1 1 1 0]; % S0+I0+R0 <= 1
beq = [1.000001];

[modeledNewInfections, modeledCumulativeDeaths, totalModel, breakthroughInfections] = simulate([ki kr kd kv kb kp S0 I0 R0 d]);

plot(modeledNewInfections, "LineWidth",line_size)
plot(modeledCumulativeDeaths, "LineWidth", line_size)
plot(breakthroughInfections, "LineWidth",line_size)
% plot(totalModel(:,4), "LineWidth",line_size) % Vaccinations, Breaks y axis legibility
legend("Mock New Infections", "Mock Cumulative Deaths", "Modeled New Infections", "Modeled Cumulative Deaths", "Modeled Breakthrough Infections")

vaxpop = totalModel(:,4);
vaxbreak = breakthroughInfections;
% Save results
save("competition.mat", "vaxpop", "vaxbreak");

function [newInfections, cumulativeDeaths, totalModel, breakthroughInfections] = simulate(input)
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
    V0 = 0;
    D0 = 1-(S0+I0+R0+V0);
    d  = round(input(10));

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
    
    % Initial state for phase 2 is end of phase 1
    initialState(1) = modeledBeforeVaccine(end,1);
    initialState(2) = modeledBeforeVaccine(end,2);
    initialState(3) = modeledBeforeVaccine(end,3);
    initialState(4) = modeledBeforeVaccine(end,4);
    initialState(5) = modeledBeforeVaccine(end,5);
    
    sys_sir_base = ss(A, zeros(5,1), eye(5), zeros(5,1),1);
    modeledAfterVaccine = lsim(sys_sir_base, zeros(d,1),linspace(0,d-1,d),initialState);
    
    totalModel = [modeledBeforeVaccine; modeledAfterVaccine];
    
    % Convert and return results
    modeledNewInfections = totalModel(:,1) * ki; % Susceptible -> infected
    breakthroughInfections = totalModel(:,4) * kb; % Vaccinated -> infected
    modeledNewVaccinations = [zeros(400-d,1); modeledAfterVaccine(:,1) * kv + modeledAfterVaccine(:,3) * kv];
    newInfections = modeledNewInfections + breakthroughInfections + modeledNewVaccinations;
    cumulativeDeaths = totalModel(:,5);

    % Assume people are infected for 2 weeks after a breakthrough infection
    window = 14; % 2 weeks
    block = ones(1,window);
    breakthroughInfections = conv(breakthroughInfections, block, "same"); % 14 running sum
end