% Trial 1 of Original SIRD Model
SIRD_1(1:4,1) = [1 0 0 0].'; % Initial State
A_1 = [0.95 0 0 0; 0.05 0.87 0 0; 0 0.12 1 0; 0 0.01 0 1]; % Original SIRD Matrix
for i = 2:200
    SIRD_1(:,i) = A_1*SIRD_1(:,i-1);
end
figure;
for j = 1:4
    plot(1:200,SIRD_1(j,:))
    hold on
end
xlabel("Time");
ylabel("x");
title("SIRD Basic Model", fontsize = 20)
legend("Susceptible","Infected","Recovered","Deceased")

% Trial 2 of Original SIRD Model
SIRD_1(1:4,1) = [1 0 0 0].'; % Initial State
A_1 = [0.9 0.1 0 0; 0.1 0.85 0 0; 0 0.1 1 0; 0 0.05 0 1]; % Original SIRD Matrix
for i = 2:200
    SIRD_1(:,i) = A_1 * SIRD_1(:,i-1);
end
figure;
for j = 1:4
    plot(1:200,SIRD_1(j,:))
    hold on
end
xlabel("Time");
ylabel("x");
title("Trial 2 of SIRD Normal Model")
legend("Susceptible","Infected","Recovered","Deceased")

%How does the output of the model tend to converge over time?
%Deceased rate converge to the ratio that I set in the matrix A.
%Susceptible, Infected, and Recovered rate converge to 0.

% SIRD Matrix when Reinfections are possible
SIRD_2(1:4,1) = [1 0 0 0].'; % Initial State
A_2 = [0.95 0 0.2 0; 
       0.05 0.89 0 0; 
       0 0.1 0.8 0; 
       0 0.01 0 1]; % Reinfections SIRD Matrix

for i = 2:1000
    SIRD_2(:,i) = A_2 * SIRD_2(:,i-1);
end
figure;
for j = 1:4
    plot(1:1000,SIRD_2(j,:))
    hold on
end
xlabel("Time");
ylabel("x");
title("SIRD when Reinfections are Possible",fontsize = 20);
legend("Susceptible","Infected","Recovered","Deceased")
% Decesed Behavior converge to 1 as time goes on.


SIRD_3(1:4,1) = [1 0 0 0].'; % Initial State
A_3 = [0.95 0.04 1 0; 
       0.05 0.85 0 0; 
       0 0.1 0 0; 
       0 0.01 0 1]; % Reinfections SIRD Matrix

for i = 2:1000
    SIRD_3(:,i) = A_3 * SIRD_3(:,i-1);
end
figure;
for j = 1:4
    plot(1:1000,SIRD_3(j,:))
    hold on
end
xlabel("Time");
ylabel("x");
title("SIRD when Reinfections are Possible")
legend("Susceptible","Infected","Recovered","Deceased")
% Decesed Behavior converge to 1 as time goes on.