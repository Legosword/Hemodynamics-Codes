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
%x = x - x0;
%P5 = P5 - x0;
disp('Solution vector x = ')
disp(x)
disp(P5)
disp(P6)
disp(['Residual after ', num2str(counts), ' iterations = ', num2str(resi)])