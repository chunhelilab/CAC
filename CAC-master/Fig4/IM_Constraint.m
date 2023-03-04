
function [C_ineq, C_eq, Cineq_grad, Ceq_grad] = IM_Constraint(x, ConstraintParams)
%x is the tunable parameters of the function. 
%This function includes all constraints on the point x, including all
%bifurcation conditions, etc. 

%Here we use this function to ensure that specified stable fixed points
%stay in existence. 

Function = ConstraintParams.Function; %This should be the main function, and it should 
%include a time component as well. 

Jacobian = ConstraintParams.Jacobian; 

original_fp = ConstraintParams.OriginalFP; 
original_params = ConstraintParams.OriginalParams;
L = ConstraintParams.L; 



 BifurcationVector = [0 0 0 1];  %A 1,0 array of which fp should be preserved.
%1: Preserve the specified state
%2. The specified state can undergo a bifurcation. 
%This vector can be modified if it is acceptable for certain states to be
%lost to bifurcation. I


%BifurcationVector = [1 1]; %%CL 2 states

% NumSS=size(original_fp,2);
% BifurcationVector = ones(NumSS,1);

l_Bifurcation = sum(BifurcationVector); 
C_ineq = zeros(1, l_Bifurcation); %%CL 
tol =ConstraintParams.Tolerance ; %The largest real eigenvalue should thus not be larger than -.05; 

inparams.in_params = original_params;
inparams.out_params = x;
 

for i=1:l_Bifurcation
    inparams.fpinit = original_fp(:, i);
    C_ineq(i) = Adaptive_FixedPointContinuer_Bifurcations(Function,Jacobian, inparams) + tol;     
end

%%CL, sparsity constraint, change S
% S=2.5;
% L=20;
% C_ineq(l_Bifurcation+1)=sum(abs(x-S))-L;
lambda=sum(abs(x-original_params));
C_ineq(l_Bifurcation+1)=lambda-L;
C_eq = []; 
Cineq_grad = []; 
Ceq_grad = []; 

end
function [maxRealEigenPart] = Adaptive_FixedPointContinuer_Bifurcations(func,J, params)
%This is a straight forward method to calculate the largest real part the
%the eigenvalues associated with a particular stable state as that state is
%continued. Since this value changes very suddenly as a bifurcation is
%appraoched, an adaptive scheme had to be implemented. 

format long; 

%What is the fixed point that we want to continue?
fpinit = params.fpinit; 
 originalfp = fpinit; 
 
%What is the initial set of parameters that corresponds to this fixed point?
in_params = params.in_params; 
num_params = size(in_params, 2); 

%What set of parameters are we going to? 
out_params=params.out_params; 

%How many steps should the continuation take? 
N1=  20; %Number of steps to do the initial continuation
N2 = 21; %Number of steps to do teh adaptive continuation of. 
 

%Calculate the line for the homotopy method. 
param_line = zeros(num_params, N1); 
for i=1:num_params
    param_line(i, :) = linspace(in_params(i), out_params(i), N1); 
end



%This N should be relativly small. 
norm_mat = zeros(16, N1+1);  %%%CL change dimension of network
norm_mat(:, 1) = fpinit;
for i=1:N1
    params_vals=param_line(:, i); 
    f = @(t, x)func(t, x, params_vals);  
    %Calculate the new fixed point
    [T, Y]=ode15s(f, [0, 2e3], fpinit);
    fpinit_test=Y(end, :)';
    norm_mat(:,i+1) = fpinit_test; %Save the change in norm between fixed point. 
    %Those places where the stable state is changing the most are where a
    %bifurcation could be occurring. 
    fpinit = fpinit_test; 
end
%This is an adaptive scheme to calculate the largest real part of the
%jacobian of the system where the stable states are changing the most. 
%This is essentially a numerical variational problem. 
%dparams = out_params - in_params; %These are the changes in each fixed point. 
dnorm = sum(abs((norm_mat(:, 2:end) - norm_mat(:, 1:end-1))));  


w = sqrt(1 + 1e10*dnorm.^2); 
alphak = [0 cumsum(w)/sum(w)];
remesh= interp1(alphak, 0:N1, linspace(0, 1, N2));
remesh(end)= remesh(end-1);
remesh = remesh./max(remesh); 

newparams = zeros(size(in_params, 2), N2); 
for i=1:size(in_params, 2)
    newparams(i, :) = spline(linspace(0, 1, N1),param_line(i, :), remesh); 
end
Eigen = -1; 
%EigMat = zeros(1, N2);
fpinit = originalfp; 
for i=1:N2
    params_vals=newparams(:, i); 
    f = @(t, x)func(t, x, params_vals);  
    jtemp = @(x)J(1, x, params_vals); 
    %Calculate the new fixed point
    [T, Y]=ode15s(f, [0, 1e3], fpinit);
    fpinit_test=Y(end, :)';
    warning off
    Eigentest = max(real(eigs(jtemp(fpinit)))); 
    Eigen = max(Eigen, Eigentest); 
    %EigMat(i) = Eigen; 
    fpinit = fpinit_test; 
end
maxRealEigenPart = Eigen;  
end
