A = [0.95 0.04 0.01 0; 0.05 0.85 0.01 0; 0 0.1 0.98 0; 0 0.01 0 1];
SIRD(1:4,1) = [1 0 0 0].';
for i = 2:200
    SIRD(:,i) = A * SIRD(:,i-1);
end
figure;
for j = 1:4
    plot(1:200,SIRD(j,:))
    hold on
end
xlabel("Time");
ylabel("x");
legend("Susceptible","Infected","Recovered","Deceased")