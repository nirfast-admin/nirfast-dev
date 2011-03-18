function [A,b,Aeq,beq,lb,ub]=constrain_lsqfit(wv,num_param)
% Assumes HbO, deoxyHb, and water are parameters in the reconstruction.


A = zeros(wv,num_param-2);
b = zeros(wv,1);

% build matrix A for Ax=b contraint 
A(1,1) = 1; % HbO
A(2,2) = 1; % deoxyHb
A(3,1:2) = 1; % Hbt
A(4,3) = 1; % water 


b(1:3,1) = 0.1; % HbO, deoxyHb, Hbt to less than 0.1 mM.
b(4,1) = 1; % water to less than 100%

for i = 1:num_param-5
   A(4+i,3+i) = 1; % exog. agents
   b(4+i,1) = 0.1; % exog. agents to less than 0.1mM
end

% set lower bounds:
lb = zeros(num_param-2,1);
ub = b(1:2,1);
ub(3,1) = 1;
for i = 1:num_param-5
   ub(3+i,1) = 0.1; % exog. agents to less than 0.1mM
end

Aeq = zeros(wv,num_param-2);
beq = zeros(wv,1);