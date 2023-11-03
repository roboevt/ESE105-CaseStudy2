
% Interaction Model without reinfections
flow1 = 0.1; % Population flow between two cities without any restrction
flow2 = 0.06; % People are less likely to travel between cities
ratio = 5; % The ratio of population x and y
SIRD(1:8,1) = [1-1/(1+ratio) 0 0 0 1/(1+ratio) 0 0 0].'; %The initial population of two cities
% A_1 is the normal infection rate(including reinfections) without virus variants
A_1 = [0.95*(1-flow1) 0                 0.3*(1-flow1)  0    0.95*flow1      0                0.3*flow1     0;       
       0.05*(1-flow1) 0.85*(1-flow1/5)  0              0    0.05*flow1      0.85*flow1/3     0             0; 
       0              0.14*(1-flow1/2)  0.7*(1-flow1)  0    0               0.14*flow1       0.7*flow1     0;
       0              0.01              0              1    0               0                0             0;
       0.95 * flow1   0                 0.3*flow1      0    0.95*(1-flow1)  0                0.3*(1-flow1) 0; 
       0.05 * flow1   0.85*flow1/5      0              0    0.05*(1-flow1)  0.85*(1-flow1/3) 0             0;   
       0              0.14*flow1/2      0.7*flow1      0    0               0.14*(1-flow1)   0.7*(1-flow1) 0;
       0              0                 0              0    0               0.01             0             1];
% Omicron has higher infection rate, lower death rate, and higher recovery
% rate
% A_2 is Omicron infection rate with no public policy
A_2 = [0.9*(1-flow1)  0                 0.3*(1-flow1)  0    0.9*flow1      0                 0.2*flow1      0;       
       0.1*(1-flow1)  0.85*(1-flow1/5)  0              0    0.1*flow1      0.85*flow1/5      0              0; 
       0              0.145*(1-flow1/2) 0.7*(1-flow1)  0    0              0.145*flow1/2     0.8*flow1      0;
       0              0.005             0              1    0              0                 0              0;
       0.9 * flow1    0                 0.3*flow1      0    0.9*(1-flow1)  0                 0.2*(1-flow1)  0; 
       0.1 * flow1    0.85*flow1/3      0              0    0.1*(1-flow1)  0.85*(1-flow1/5)  0              0;   
       0              0.145*flow1       0.7*flow1      0    0              0.145*(1-flow1/2) 0.8*(1-flow1)  0;
       0              0                 0              0    0              0.005             0              1];

% A_3 is Omicron infection rate with travel policy (people can't travel if
% you are infected)
A_3 = [0.9*(1-flow1)  0                 0.3*(1-flow1)  0    0.9*flow1      0                 0.2*flow1      0;       
       0.1*(1-flow1)  0.85              0              0    0.1*flow1      0                 0              0; 
       0              0.145*(1-flow1/2) 0.7*(1-flow1)  0    0              0.145*flow1/2     0.8*flow1      0;
       0              0.005             0              1    0              0                 0              0;
       0.9 * flow1    0                 0.3*flow1      0    0.9*(1-flow1)  0                 0.2*(1-flow1)  0; 
       0.1 * flow1    0                 0              0    0.1*(1-flow1)  0.85              0              0;   
       0              0.145*flow1       0.7*flow1      0    0              0.145*(1-flow1/2) 0.8*(1-flow1)  0;
       0              0                 0              0    0              0.005             0              1];

% A_4 is normal infection rate with naturally low population flow
A_4 = [0.95*(1-flow2) 0                 0.3*(1-flow2)  0    0.95*flow2      0                0.3*flow2     0;       
       0.05*(1-flow2) 0.85*(1-flow2/5)  0              0    0.05*flow2      0.85*flow2/3     0             0; 
       0              0.14*(1-flow2/2)  0.7*(1-flow2)  0    0               0.14*flow2       0.7*flow2     0;
       0              0.01              0              1    0               0                0             0;
       0.95 * flow2   0                 0.3*flow2      0    0.95*(1-flow2)  0                0.3*(1-flow2) 0; 
       0.05 * flow2   0.85*flow2/5      0              0    0.05*(1-flow2)  0.85*(1-flow2/3) 0             0;   
       0              0.10*flow2/2      0.7*flow2      0    0               0.14*(1-flow2)   0.7*(1-flow2) 0;
       0              0                 0              0    0               0.01             0             1];
% Omicron has lower infection rate, higher death rate, and lower recovery
% rate

% A_5 is Delta infection rate without public policy
A_5 = [0.98*(1-flow2)  0                 0.3*(1-flow2)  0    0.98*flow2      0                 0.3*flow2     0;       
       0.02*(1-flow2)  0.88*(1-flow2/6)  0              0    0.02*flow2      0.88*flow2/6      0             0; 
       0               0.105*(1-flow2/2) 0.7*(1-flow2)  0    0               0.105*flow2/2     0.7*flow2     0;
       0               0.015             0              1    0               0                 0             0;
       0.98 * flow2    0                 0.3*flow2      0    0.98*(1-flow2)  0                 0.3*(1-flow2) 0; 
       0.02 * flow2    0.88*flow2/6      0              0    0.02*(1-flow2)  0.88*(1-flow2/6)  0             0;   
       0               0.105*flow2/2     0.7*flow2      0    0               0.105*(1-flow2/2) 0.7*(1-flow2) 0;
       0               0                 0              0    0               0.015             0             1];

% A_5 is Delta infection rate with public policy that susceptible and
% infected people can't travel
A_6 = [0.98  0                 0.3            0    0         0                 0.3           0;       
       0.02  0.88              0              0    0         0                 0             0; 
       0     0.105*(1-flow2/2) 0.7*(1-flow2)  0    0         0.105*flow2/2     0.7*flow2     0;
       0     0.015             0              1    0         0                 0             0;
       0     0                 0              0    0.98      0                 0             0; 
       0     0                 0              0    0.02      0.88              0             0;   
       0     0.105*flow2/2     0.7*flow2      0    0         0.105*(1-flow2/2) 0.7*(1-flow2) 0;
       0     0                 0              0    0         0.015             0             1];

for i = 2:200  
    SIRD(:,i) = A_1 * SIRD(:,i-1);
end
for i = 201:300
    SIRD(:,i) = A_2 * SIRD(:,i-1);
end
for i = 301:400
    SIRD(:,i) = A_3 * SIRD(:,i-1);
end
for i = 401:600
    SIRD(:,i) = A_4 * SIRD(:,i-1);
end
for i = 601:700
    SIRD(:,i) = A_5 * SIRD(:,i-1);
end
for i = 701:800
    SIRD(:,i) = A_6 * SIRD(:,i-1);
end
% Compare rate per population
SIRD(1:4,:)/ratio;
figure;
for j = 1:8
    plot(1:800,SIRD(j,:))
    hold on
end
xlabel("Time");
ylabel("x");
title("Travel Model in Different Phases",fontsize = 20)
legend("Susceptible 1","Infected 1","Recovered 1","Deceased 1","Susceptible 2","Infected 2","Recovered 2","Deceased 2")

SIRD_susceptibleDifference = SIRD(1,:) - SIRD(5,:);
SIRD_infectedDifference = SIRD(2,:) - SIRD(6,:);
SIRD_recoveredDifference = SIRD(3,:) - SIRD(7,:);
SIRD_deceasedDifference = SIRD(4,:) - SIRD(8,:);

figure;
tiledlayout(2,2)
% Topleft plot
nexttile
plot(1:800,SIRD_susceptibleDifference);
title('Difference of Susceptible People')

% Topright plot
nexttile
plot(1:800,SIRD_infectedDifference);
title('Difference of Infected People')

nexttile
plot(1:800,SIRD_recoveredDifference);
title('Difference of Recovered People')

nexttile
plot(1:800,SIRD_deceasedDifference);
title('Difference of Deceased People')

