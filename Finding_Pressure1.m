%------------------------------%
%This program requires package tool box.
%Goal of script: Uses Power law model to calculate Pressures for 6 nodes
% Pressure 6 is labeled as P6, is our initial guess. Q0 is another choice,
% that is needed to solve Pressures at other nodes.
% Uses Multidimentional Newton's method to solve for P_i for i = 1,...,4
% Initial guess is vector x.
%Run Finding_Pressure1.m, then run Finding_Pressure.m
%------------------------------%
syms P1 P2 P3 P4
%Grab Radius r from data
data = readtable('testData.csv');
R = data.Var5;
format long 
%some Q0 to start it off and solve
Q0= 1;
%P = [P1 P2 P3 P4 P5 P6];
%our guess to solve the system
P6 = 0;
n = 1;
%length L
L = 1;
k = 3;
%coeffient of R in Q
Constant = (pi*n)/(3*n+1)*(1/(2*L*k))^(1/n);
M = zeros(1,length(R));
%disp(M)
for i = 1:length(R)
    %Q at every edge without Pressure
   M(i) = Constant*(R(i))^((3*n+1)/(n));
end
P5 = (8*L*k)/(pi*R(7)^4);
%disp ("This is P5 =" + P5)

% count, epsilon, maxcount,res
eps = 1e-6;
countMax = 1000;
count = 0;
res = 10;  % Initial residual to start the loop
% initial guesses for [P1, P2, P3, P4]
x = [5*P5; 4*P5; 2*P5; 3*P5];
%x0 = mean(x);
%x0 = 0;
%x=x+x0;
%P5 = P5 + x0;
% function handles for each equation
F = cell(4,1);
F{1} = @(P1, P2) M(1)*(P1 - P2)^(1/n) - Q0;
F{2} = @(P2, P3, P4) M(2)*(P2 - P4)^(1/n) + M(3)*(P2 - P3)^(1/n) - Q0;
F{3} = @(P2, P3, P4) M(3)*(P2 - P3)^(1/n) + M(4)*(P4 - P3)^(1/n) + M(5)*(P4 - P5)^(1/n) - Q0;
F{4} = @(P3, P4) M(5)*(P4 - P5)^(1/n) + M(6)*(P3 - P5)^(1/n) - Q0;
% Jacobian matrix
J = jacobian([F{1}(P1, P2), F{2}(P2, P3, P4), F{3}(P2, P3, P4), F{4}(P3, P4)], [P1, P2, P3, P4]);
% Convert the Jacobian and functions to function handles
F_numeric = @(x) [F{1}(x(1), x(2)); F{2}(x(2), x(3), x(4)); F{3}(x(2), x(3), x(4)); F{4}(x(3), x(4))];
J_numeric = matlabFunction(J, 'Vars', {P1, P2, P3, P4});
% multi-dimentional newton's method
while res > eps && count < countMax
    % Evaluate F and J at the current guess x
    F_val = F_numeric(x);
    J_val = J_numeric(x(1), x(2), x(3), x(4));
    % solve for update vector delta x
    delta_x = -J_val \ F_val;
    
    % Update x
    x = x + delta_x;
    % update residual and count
    res = norm(F_val);
    count = count + 1;
    % couldddddd display, not really needed.
    fprintf('Iteration %d: Residual = %e\n', count, res);
end
% Display solution
%x = x - x0;
%P5 = P5 - x0;
disp('Solution vector x = ')
disp(x)
disp(P5)
disp(P6)
disp(['Residual after ', num2str(count), ' iterations = ', num2str(res)])