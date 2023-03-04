%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is an implementation of the adaptive minimum action method developed
%by Zhou et al. in Zhou X, Ren W, E W (2008) An adaptive minimum action method for the study
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [fval, phi] = AMAM(params, funcs, dfuncs)
format long

%The Jacobian should be such that it can take in a 2-d array of values
%through time, and return the 3-d time Jacobian. Time should be on the
%third axis. 


xffp = params.init;
%The fixed points should be given in columns. 
tmax = params.TMax;
n = params.N; 
%Note that n is the number of points to use, not the dimension.
m = params.Dimension;
c = params.c;
maxiter = params.MaxIter;
kmax = params.K;
remesh = params.Remesh;
qmin = params.q;
phi = params.PhiInit;

phil = xffp(:,1); 
phir = xffp(:,2);

T = -tmax/2:tmax/n:tmax/2;
delt = tmax/n*ones(m, n);       
k = 0;  

%If an initial path is not set, use a straight line. 
if isempty(params.PhiInit) 
    phi = pathfinder(phil, phir, T);
end

while k < kmax 

    %ub is alway set to be inf. 
     ub = inf*ones(m,n-1);
     lb = -inf*ones(m, n-1); 
      

     if params.ObjGrad==1     
         options = optimset('Algorithm', 'interior-point','Hessian',{'lbfgs',5},...
        'MaxIter',maxiter,  'MaxFunEvals', 3000000,...
    'AlwaysHonorConstraints', 'none','Display', 'off', 'GradObj','on');
        func = @(x)Sg(x, delt,funcs,dfuncs,phir, phil);
        
     else
         options = optimset('Algorithm', 'interior-point','Hessian',{'lbfgs',5},...
        'MaxIter',maxiter,  'MaxFunEvals', 3000000,...
       ...% 'TolX', 1e-10,'TolCon', 1e-10, 'TolFun', 1e-10,
    'AlwaysHonorConstraints', 'none','Display', 'off');
        func = @(x)S(x, delt,funcs,phir, phil); 
       
     end

    phi = phi(:, 2:end-1);



    [x2, fval] = fmincon(func, phi, [],[],[],[],...
                         lb,ub,[],options);

     x2mod = [phil, x2, phir];
     w = monitor(x2mod, delt, c);
     q = qfunc(w, delt);
     if remesh == 1 && q > qmin
         [T, phi, delt] = gridremesh(n,m,w,delt,T, x2mod);
         else
         phi = x2mod;
     end
     k = k+1  ; 
end
phi = phi(:, 2:end-1); 
         
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [out, g] = Sg(phi, delt, func,dfunc, phirp, philp)
%Approximate the action functional according to the midpoint
%rule. The time derivative of phi is taken inside the function.

%Create the relevant phi vectors
phi = [philp,phi, phirp];
phir = phi(:,2:end);
phil = phi(:,1:end-1);
phihalf = (phir + phil)/2;




%Perform the trapezoidal approximation. 
summand = sum(((phir - phil)./(delt) - func(phihalf)).^2);
%summand = sum(summand.^2);             
out = .5*sum(summand.*delt(1, :));

g = gradS(phi(:,2:end-1) , delt, func, dfunc, phirp, philp); 
end

function  out = S(phi, delt, func, phirp, philp)
%Approximate the action functional according to the midpoint
%rule. The time derivative of phi is taken inside the function.

%Create the relevant phi vectors
phi = [philp,phi, phirp];
phir = phi(:,2:end);
phil = phi(:,1:end-1);
phihalf = (phir + phil)/2;




%Perform the trapezoidal approximation. 
summand = sum(((phir - phil)./(delt) - func(phihalf)).^2);
%summand = sum(summand.^2);             
out = .5*sum(summand.*delt(1, :));



end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = monitor(phi, delt, c)
%The monitor function
format long
phir = phi(:,2:end);
phil = phi(:,1:end-1);
phit = (phir - phil)./delt;
derivmagnitude = sum(phit.^2);
w = sqrt(1 + c*derivmagnitude);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Tnew, phi, delt] = gridremesh(n,m,w,delt,T,x2mod)
%The adaptive remeshing procedure, as in E et al. 
format short
%Calculated the integral of w using the midpoint rule. 
intw = sum(w.*delt(1, :));
%Calculate alpha in accordance with the paper.
alphak = [0,cumsum((w.*delt(1, :))/intw)];
%The uniform alpha. 
alphaold = 0:1/n:1;
%This is a linear interpolation of T as a function of alpha from the
%stretched alphak to the uniform alpha. 
Tnew = interp1(alphak, T, alphaold, 'linear');
%interp1 was giving a NaN in this position for reasons unknown, since it
%should just give the end value back. This is an ragged solution to the
%problem. 
Tnew(end) = -T(1);
%The matrix of delt. 
delt = ones(m, 1)*(Tnew(2:end) - Tnew(1:end-1));
%Find the new phi by a cubic interpolation over the remeshed T.
phi = spline(T, x2mod, Tnew);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = qfunc(w, delt)
%The calculation of q, in accordance with E et al. 
       maxx = max(w.*delt(1, :)); minx = min(w.*delt(1, :));
       q = maxx/minx;
       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function path = pathfinder(start, ending, time)
%given to starting points, creates a matrix of a linear transit between
%them. 

dims = size(start);
dim = dims(1);
time = time + time(end);

% [xstart, ystart] = meshgrid(start);
% [xend, yend] = mehsgrid(end);
% [xt, yt]= meshgrid(time);
% 
% path

xstart = start*ones(1, length(time));
diff = (ending- start)*ones(1, length(time));
timep = (time/(max(time)))';
times = (timep*ones(1, dim))';
path = xstart + diff.*times;


end

function out = gradS(phi, delt, func, dfunc, phirp, philp)

%Add down, not across, for the gradient. 
phi = [philp,phi, phirp];
phir = phi(:,2:end);
phil = phi(:,1:end-1);
phihalf = (phir + phil)/2;
%These go across all dimensions and time.
sz = size(phihalf); 
% out = zeros(sz(1), sz(2)-1); 
% ndim = sz(1); 
%For now, we will have to use a single loop
T0 = ((phir - phil)./(delt) - func(phihalf)); 

T0p = reshape(T0, [1,sz(1), sz(2)]); 
deltp = reshape(delt,  [1, sz(1),  sz(2)]); 
phijacob = dfunc(phihalf);

out = zeros(1, sz(1), sz(2)-1); 
for i=1:sz(1) 
    a1 = T0p(1, i, 1:end-1); 
    a2 = T0p(1, i, 2:end); 
    
    step1 = a1.*(1./(deltp(1, i, 1:end-1)) - 1/2*phijacob(i,i,1:end-1)).*deltp(1, i, 1:end-1);
    step2 = a2.*(-1./(deltp(1, i, 2:end)) - 1/2*phijacob(i,i,2:end)).*deltp(1, i, 2:end);
    vals = 1:sz(1); 
    vals(i) = []; 
    step3 = -1/2*(sum( T0p(1, vals, 1:end-1).*permute(phijacob(vals,i, 1:end-1), [2, 1, 3]).*deltp(1, vals, 1:end-1), 2) + sum( T0p(1, vals, 2:end).*permute(phijacob(vals,i, 2:end), [2, 1, 3]).*deltp(1, vals, 2:end), 2)); 

    out(1, i, :) = step1 + step2 + step3; 
    
end

out = squeeze(out); 


end


function AMAM_Params=AMAM_params(Parameters)
%These parameter can be modified if desired, but they are not included for
%modification above. 
AMAM_Params={}; 
AMAM_Params.N=100; %The number of points in the minimum action path. 
%Larger is more accurate, but slower. 
AMAM_Params.TMax = 1; %The time range of the minimum action path.%20 
%Larger values can be more accurate, but can also lead to instabilities. 
AMAM_Params.Dimension = size(Parameters.StableStates, 1); %The dimension of the system. 
AMAM_Params.c = 1e10; %The remeshing parameter, c. Larger values give greater remeshing. 
AMAM_Params.MaxIter = 300; %The number of optimization steps to run between remeshing steps. 
AMAM_Params.K= 2; %The number of total remeshing steps to do; 
AMAM_Params.Remesh = 1; %Turn the remesh on or off. If off, the algorithm will be default perform K*MaxIter itereations. 
AMAM_Params.q = 3; %The q parameter from the original paper. Is a measure of path smoothness. 
AMAM_Params.ub = []; %If desired, constraints can be set of the path as well. This is not done in our manuscript. 
AMAM_Params.lb = []; %As above, lower bounds can be set on the path as well. This is also not done at all in our manuscript. 
AMAM_Params.PhiInit =[]; %The default intial path between two states is a straight line. This can be modified here. 
AMAM_Params.ObjGrad = Parameters.JacobianPresent;%We require a Jacobian in this implementation. 
    
end