% The following matrix implements the SIR dynamics example from Chapter 9.3
% of the textbook.
A = [0.95 0.04 0 0; 0.05 0.85 0 0; 0 0.1 1 0; 0 0.01 0 1];

% The following matrix is needed to use the lsim function to simulate the
% system in question
B = zeros(4,1);

% initial conditions (i.e., values of S, I, R, D at t=0).
x0 = [0.9 0.1 0 0];

% Here is a compact way to simulate a linear dynamical system.
% Type 'help ss', 'help lsim', etc., to learn about how these functions work!!
sys_sir_base = ss(A,B,eye(4),zeros(4,1),1);
Y = lsim(sys_sir_base,zeros(1000,1),linspace(0,999,1000),x0);

% plot the output trajectory
plot(Y);
legend('S','I','R','D');
xlabel('Time')
ylabel('Percentage Population');