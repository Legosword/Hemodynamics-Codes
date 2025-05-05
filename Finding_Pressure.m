%------------------------------%
%Run this Program second
%------------------------------%
n = 0.7;
%coeffient of R in Q
Constant = (pi*n)/(3*n+1)*(1/(2*L*k))^(1/n);
M = zeros(1,length(R));
%disp(M)
for i = 1:length(R)
    %Q at every edge without Pressure
   M(i) = Constant*(R(i))^((3*n+1)/(n));
end
% count, epsilon, maxcount,res
epsi = 1e-6;
countMaxi = 100;
counts = 0;
resi = 100;  % Initial residual to start the loop
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
while resi > epsi && counts < countMaxi
    % Evaluate F and J at the current guess x
    F_val = F_numeric(x);
    J_val = J_numeric(x(1), x(2), x(3), x(4));
    % solve for update vector delta x
    delta_x = -J_val \ F_val;
    
    % Update x
    x = x + delta_x;
    % update residual and count
    resi = norm(F_val);
    counts = counts + 1;    
    % couldddddd display, not really needed.
    fprintf('Iteration %d: Residual = %e\n', counts, resi);
end

% Display solution
disp('Solution to pressures are ')
fprintf('P1 = %.5f\n', x(1));
fprintf('P2 = %.5f\n', x(2));
fprintf('P3 = %.5f\n', x(3));
fprintf('P4 = %.5f\n', x(4));
fprintf('P5 = %.5f\n', P5);
fprintf('P6 = %.5f\n', P6);

%This section is to find the Q's of the example capillary network

%Q1 = pi*n/(3*n+1)*R(1,1)^((3*n+1)/n)*(P1-P2)/(2*L*K)^(1/n);
%Q2 = pi*n/(3*n+1)*R(1,2)^((3*n+1)/n)*(P1-P2)/(2*L*K)^(1/n);
%Calculate Q with Power-Law Model.
Q1 = M(1).*(x(1)-x(2)).^(1/n);
Q2 = M(2).*(x(2)-x(4)).^(1/n);
Q3 = M(3).*(x(2)-x(3)).^(1/n);
Q4 = M(4).*(x(4)-x(3)).^(1/n);
Q5 = M(5).*(x(4)-P5).^(1/n);
Q6 = M(6).*(x(3)-P5).^(1/n);
Q7 = M(7).*(P5);
disp('Solutions to Volumetric flow rates')
fprintf('Q1 = %.5f\n', Q1);
fprintf('Q2 = %.5f\n', Q2);
fprintf('Q3 = %.5f\n', Q3);
fprintf('Q4 = %.5f\n', Q4);
fprintf('Q5 = %.5f\n', Q5);
fprintf('Q6 = %.5f\n', Q6);
fprintf('Q7 = %.5f\n', Q7);

%disp(['Residual after ', num2str(counts), ' iterations = ',
%num2str(resi)]) 55