%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%addpath('/home1/zhchen/program/czhprogram/OLAC1')
%Define the parameter struct on which the algorithm will operate. 
Parameters = {}; 
%1. Define the (controllable dynamical system). 
%This should be a function that takes in three parameters, in this order: 
%(t, x, P), where t is the time (for non-autonomous dynamical systems), x
%is the state vector, and P is the parameter vector for the system. 
Function  = @(t, x, P) Force5i(t, x, P);
Parameters.Function = Function;
%2. The states of the system must be previously identified to run OLAC.
N=20;
load stable_4ss.mat; %%default
load y_4ss.mat; %%default
%load SSa03b08sa04.mat; %%default
%load good_param/stable_7.mat; %%default
%load good_param/y_7.mat; %%default

StableStates=stable;
Parameters.StableStates = StableStates;

%3. OLAC identifies the optimal parameter combination for a given goal,
%which is a function that evaluates the action of the dynamical system.
%That function must be defined for OLAC. 
%This function should be one which takes in the action for each transition
%in a dynamical system, and returns a single value. 
%For example, in the paper we use the occupancy of a single state for this
%value, as well as particular transition rates.  
OLAC_Functional = @(ActionArray, StableStatesInExistence)IM_OLACFunctional(ActionArray, StableStatesInExistence); 
Parameters.Functional = OLAC_Functional;  


%4. Stable State Identifier
%In the course of the simulation bifurcations can occur and two states can
%become one. A function is required to distinguish different stable 
%states. 
ID_Function = @(x)MISE_SSC(x, StableStates); 
Parameters.ID_Function = ID_Function; 
%5. Input the (controllable) Jacobian of the dynamical system. 
%This should take in the time, (t), state vector (x), and the set of
%parameters, and return the square matrix that is the Jacobian of the
%dynamical system. 

%Importantly, the Jacobian MUST be able to be vectorized - it must be able
%to be evaluated over a time series of different states of the system. 
%See the included example for details. 
%While this Jacobian is not technically required to run OLAC, the speed up 
%when it is present is immense (>100x), making it effectively required for 
%any realistic examples. 

Jacobian = @(t, x, P) IM_Jacobian(t, x, P); 
Parameters.Jacobian = Jacobian;
Parameters.JacobianPresent = 1; %Set this to 0 if no analytic Jacobian is given. 
NumP=189;%%number of parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OLAC Parameters
%Parameters to tune the main optimization step in OLAC.
%These parameters might need tuning for optimal performance. 
OLAC_Parameters = {}; 

OLAC_Parameters.MaxIter = 3; %Set the maximum number of iterations for Olac. 
OLAC_Parameters.ObjectiveLimit = -2; %Set a tolerance for the objective function, defined above. 
OLAC_Parameters.ub = []; %Bound parameter can be included to to maintain physiologically realistic parameter values. 
OLAC_Parameters.lb=0*ones(1,NumP); %lower bound for the parameters. 
%These bounds can be modified to achieve better convergence.

% OLAC_Parameters.InitParams =2.5*ones(1,NumP); %The initial parameter set for the model. 
a=10;n=2;d=0.3;l=800;m=10;u=0.2;
%{
y0=[a	a	 a	a	d	d	d  d	u	 24,...      %1-10
    m	m	m m	m	m	m	4	n	4,...     %11-20
    n	n    n	n	n	n	l	l	 1000 	l,...  %21-30
    l l l  l l  4.3 l 0.02 1.5  m,... %31-40
    3000 4 2 2000 n l l l l l,...                 %41-50
    a a a a a a a a a a,...                           %51-60
    a a d d d d d d d d,...                             %61-70
    d d d d m m m m u u,...                      %71-80
    15 m m m u m m m m m,...                  %81-90
    15 m m m m m m m u m,...                  %91-100
    15 m u u m m u m 120 15,...               %101-110                 
    m m n n n n n n n n,...                           %111-120
    n n n n n n n n n n,...                             %121-130
    n n n n 4 n n n n n,...                             %131-140
    n n n n n n n n n n,...                             %141-150800
    l,l,l,l,l,l,l,l,l,l,...         %151-160
    l,l,l,l,l,l,l,l,l,l,...         %161-170
    l,l,l,l,20,l,l,l,l,l,...         %171-180
    l,l,l,l,u,n,l,m,n];                                       %181-190
%}
%load y_4ss.mat
L=0.6;%%%intervention strength L
OLAC_Parameters.InitParams =y;
Parameters.OLAC_Parameters = OLAC_Parameters;
%Make iterative figures? Note that this function is not modular and will
%need to be adapted to a different model, if desired. 
Parameters.MakeFigures = 0; 
%Parameters.MakeFigures = 0; 
%%%%%%%%%%%%%%%%%%%%%%%%%Bifurcation Constraint Function%%%%%%%%%%%%%%
%Bifurcation Constraint Function
%6.
%A nonlinear constraint should be given that gives conditions about how
%close each stable state can get to bifurcation. 
%Other constraints of both inequality and equality type, can be stipulated
%as well, in accordance with the documentation of fmincon. 
NonlinearConstraintParams.L = L;
NonlinearConstraintParams.Function = Function;
NonlinearConstraintParams.Jacobian = Jacobian; 
NonlinearConstraintParams.OriginalFP = StableStates;
NonlinearConstraintParams.OriginalParams = OLAC_Parameters.InitParams;
NonlinearConstraintParams.Tolerance = 0.4; %This value can be modified 
%to ensure proper functioning of the eigenvalue check procedure. 
ConstraintFunction = @(x)IM_Constraint(x, NonlinearConstraintParams); %CL
%This mixing time constraint used in the example. 
%ConstraintFunction = @(x)MISE_Constraint_MixingTime(x, Parameters); 
Parameters.ConstraintFunction = ConstraintFunction; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now execute OLAC. 
%With all of the parameters set up, the execution is easy. 
%EXECUTE!
OutPut = OLACFunc(Parameters); 
%OutParameters is a struct which includes the following fields: 
OutPut.ObjVal %What is the ultimate objective value achieved with
%OLAC. 
OutPut.Params(:) %What is the set of parameter values that achieves
OutPut.History.ObjectiveValues(:)
%OutPut.History.Params{:}%Give the full history of the optimization/
save('OutPut_4ss','OutPut')
