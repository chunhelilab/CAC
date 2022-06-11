function [xx,sigma,n,ycell,action]=Solver(cycle_index,par,d)
%%cycle_index: the number of random initial conditions to the ODEs to be solved
%%par: the parameters of the ODE
%%signal: Signal parameters
%%d: the diffusion coefficient 
%%%parameters for calculating paths
params={};
params.N=100; %The number of points in the minimum action path.
params.TMax = 10; %The time range of the minimum action path.
%Larger values can be more accurate, but can also lead to instabilities.
params.Dimension = 16; %The dimension of the system.
params.c = 1e14; %The remeshing parameter, c. Larger values give greater remeshing.
params.MaxIter = 300; %The number of optimization steps to run between remeshing steps.
params.K= 2; %The number of total remeshing steps to do;
params.Remesh = 1; %Turn the remesh on or off. If off, the algorithm will be default perform K*MaxIter itereations.
params.q = 3; %The q parameter from the original paconstraints can be set of the path as well. This is not done in our per. Is a measure of path smoothness.
params.ub = []; %If desired, manuscript.
params.lb = []; %As above, lower bounds can be set on the path as well. This is also not done at all in our manuscript.
params.PhiInit =[]; %The default intial path between two states is a straight line. This can be modified here.
%params.ObjGrad = Parameters.JacobianPresent;%We require a Jacobian in this implementation.
params.ObjGrad=1;
Num=16;
  %%the dimension od the system
xx=zeros(cycle_index,Num);

%%Solve odes from different initial values
for i=1:cycle_index
    i
    x0=3000*rand(1,Num);
    [t,x]=ode45(@(t,x)Force5i(t,x,par),[0,1000],x0);
    newx=x(end,:);
    x=inf*ones(1,Num);
    while norm( x(end,:)-newx(end,:) ,2 )>1e-7
        x=newx;
        [t,newx]=ode45(@(t,x)Force5i(t,x,par),[0,1],x(end,:));
    end
    xx(i,:)=newx(end,:);
end

%%Finding the stable points
 for q=1:(cycle_index-1)
     for p=(q+1):cycle_index
         if norm(xx(q,:)-xx(p,:),'fro')<10^-3
             xx(p,:)=xx(q,:);
         end
     end
 end
stable_point=unique(xx(:,:),'rows');
n=zeros(1,2);
sigma=zeros(size(xx,1),size(xx,2)^2);
for i=1:size(stable_point,1)
    [m]=find(xx(:,1)==stable_point(i,1));
    if length(m)>=1
        disp(strcat(num2str(stable_point(i,:)),' repeat ',num2str(length(m)),' times',' the location in the row xx is' ,mat2str(m)))
    end
    n(i,1)=m(1);
    n(i,2)=length(m);
    %%%calculate the covariance of each stable state
    sig=calculate_sigma(1,xx(m(1),:)',par,Num,d)';  
    for j=1:length(m)
        sigma(m(j),:)=sig;
    end
end

%Arrange the index n from large to small by fourth elements
m=size(n,1);
if(m~=1)
    for i=1:m
        tran=n(i,:);
        flag=i;
        if( i~=m )
            for j=i+1:m
                if( xx(n(j,1),4)>xx(tran(1),4) )
                    tran=n(j,:);
                    flag=j;
                end
            end
        end
        n(flag,:)=n(i,:);
        n(i,:)=tran;
    end
end
%%stable state
SS=zeros(size(xx,2),size(n,1));
for i=1:size(n,1)
    SS(:,i)=xx(n(i,1),:)';
end

Func=@(x)Force5i(1,x,par);  %Force
dFunc=@(x)jacobian(1,x,par); %Jacobian
feasi_n=size(n,1);
action=zeros(feasi_n,feasi_n);
ycell=cell(feasi_n,feasi_n);
for i=1:feasi_n
    for j=1:feasi_n
        if i==j
            action(i,j) = Inf;
        else
            Spot1 = i
            Spot2 = j%initial and end point
            params.init = SS(:, [Spot1,Spot2]);
            [ActionVal, Path] = AMAM(params, Func, dFunc);
            action(i,j) = ActionVal;
            ycell(i,j) = {[SS(:,i) Path SS(:,j)]};
        end
%         (i-1)*feasi_n+j
    end
end
end


%AMAM.m
%This is an implementation of the adaptive minimum action method developed
%by Zhou et al. in 
%Zhou X, Ren W, E W (2008) An adaptive minimum action method for the study
%of rare events. J. Chem Phys 128:104111. 
%It is to accompany the software suite with top level
%"OptimalLeastActionControl_Top.". All license information in that file
%pertains to this file as well. 
%This is sparsely documented. A more thorough description of the method can
%be found in the above paper. 
%For use with Optimal Least Action Control, nothing in this code should be
%altered or changed, unless a bug is found. 

%Author: Daniel K. Wells 
%Date: 11/11/2014
%Ver 1.0
%Email: dannykwells@gmail.com
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
delt = tmax/n*ones(m, n);%10^-3
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
       maxx = max(w.*delt(1, :)) 
       minx = min(w.*delt(1, :))
       q = maxx/minx
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
AMAM_Params.TMax = 0.01; %The time range of the minimum action path. 
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

function makefigure


load('OutPut.mat');
Parameters=OutPut.History.Parameters;
FP = Parameters.StableStates; 
out=OutPut.History.Params{2};
Params = Parameters.OLAC_Parameters.InitParams; 
Function = Parameters.Function; 
Jacobian = Parameters.Jacobian;
Func = @(x)Function(1, x, out); 
dFunc = @(x)Jacobian(1, x, out); 
szFP= size(FP, 2); 
%modification above. 
NumP=189;
AMAM_Params={}; 
AMAM_Params.N=100; %The number of points in the minimum action path. 
%Larger is more accurate, but slower. 
AMAM_Params.TMax = 0.01; %The time range of the minimum action path. 
%Larger values can be more accurate, but can also lead to instabilities. 
AMAM_Params.Dimension = size(Parameters.StableStates, 1); %The dimension of the system. 
AMAM_Params.c = 1e10; %The remeshing parameter, c. Larger values give greater remeshing. 
AMAM_Params.MaxIter = 10; %The number of optimization steps to run between remeshing steps. 
AMAM_Params.K= 1; %The number of total remeshing steps to do; 
AMAM_Params.Remesh = 3; %Turn the remesh on or off. If off, the algorithm will be default perform K*MaxIter itereations. 
AMAM_Params.q = 3; %The q parameter from the original paper. Is a measure of path smoothness. 
AMAM_Params.ub = []; %If desired, constraints can be set of the path as well. This is not done in our manuscript. 
AMAM_Params.lb =zeros(1,NumP); %As above, lower bounds can be set on the path as well. This is also not done at all in our manuscript. 
AMAM_Params.PhiInit =[]; %The default intial path between two states is a straight line. This can be modified here. 
AMAM_Params.ObjGrad = Parameters.JacobianPresent;%We require a Jacobian in this implementation. 
colors = {[203/255, 50/255, 219/255] ,... 
    [35/255, 222/255, 52/255], [26/255, 157/255, 201/255],...
    [227/255 124/255 47/255]};
place = 1; 
StorageCell = cell(2, szFP^2-szFP); 
ActionMat  =zeros(szFP); 
FP_Cell = cell(1, szFP);
close all
for ii=1:size(FP, 2)
    for jd=1:size(FP, 2)
        if ii~=jd
           
        ContinuerParams.fpinit = FP(:, ii); 
        ContinuerParams.in_params = Params; 
        ContinuerParams.out_params = out;
        FPii = FixedPointContinuer(Function, ContinuerParams);        
       
        ContinuerParams.fpinit = FP(:, jd); 
        FPjj = FixedPointContinuer(Function, ContinuerParams);     
        FP_Cell{ii} = FPii; 
        AMAM_Params.init = [FPii, FPjj]; 
        AMAM_Params.TMax = 0.01;
        AMAM_Params.N = 100; 
        AMAM_Params.Remesh = 3; 
        AMAM_Params.K= 1;
        [ActionVal, Path] = AMAM(AMAM_Params, Func, dFunc);
        StorageCell{1,place} = Path; 
        StorageCell{2, place} = ActionVal; 
        StorageCell{3, place}=[ii,jd];
        place = place +1;
        ActionMat(ii,jd) = ActionVal; 
        end

    end   
end

Num = place-1;
NewStorageCell = StorageCell; 
Places2Delete = [];
for ii=1:Num
    for jj=1:Num          
        if (ActionMat(StorageCell{3, ii}(1),StorageCell{3,ii}(2) ) + ActionMat(StorageCell{3, jj}(1),StorageCell{3,jj}(2) )) < ActionMat(StorageCell{3, ii}(1),StorageCell{3,jj}(2) )
            for kk=1:Num
                if StorageCell{3, kk} == [StorageCell{3, ii}(1) StorageCell{3, jj}(2)]
                    Places2Delete = [kk Places2Delete];  
                end
            end
        end
    end
end
for ii=1:Num
    if StorageCell{2, ii} == 0
        Places2Delete = [ii Places2Delete]; 
    end
end

NewStorageCell(:, unique(Places2Delete))=[];
Num = size(NewStorageCell, 2); 

ZEB=0:500:6000;
oct4=0:40:480;
Mdm2=0:15:180;
miR145=0:200:2400;
miR200=0:200:2400;
miR34=0:200:2400;
p53=0:200:2400;
rkip=0:200:2400;
let=0:200:2400;
lin=0:200:2400;
bach=0:200:2400;
MEK=0:200:2400;
ERK=0:200:2400;
CEBP=0:200:2400;

step(4)=50;
step(16)=15;
snail=0:step(4):step(4)*12; 
PPAR=0:step(16):step(16)*12;
[ze,oc] = meshgrid(ZEB,oct4);
[md,sn] = meshgrid(Mdm2,snail);
[mi1,mi2] = meshgrid(miR145,miR200);
[mi3,p5] = meshgrid(miR34,p53);
[rk,le] = meshgrid(rkip,let);
[li,ba] = meshgrid(lin,bach);
[me,er] = meshgrid(MEK,ERK);
[ce,pp] = meshgrid(CEBP,PPAR);
cf = @(x,y)Func([ze(:)';oc(:)';md(:)';x;mi1(:)';mi2(:)';mi3(:)';p5(:)';rk(:)';le(:)';li(:)';ba(:)';me(:)';er(:)';ce(:)';y;]);
arrows = cf(sn(:)', pp(:)');

figure(1)
h(1) = quiver(sn,pp,reshape(arrows(4, :), 13, 13), reshape(arrows(16, :),13,13), 'LineWidth', 1.2,...
    'Color',[72/255, 92/255, 99/255] ); hold on;
axis([0 step(4)*12+100 0 step(16)*12+100])

for ii=1:Num
    plot(NewStorageCell{1,ii}(4, :), NewStorageCell{1,ii}(16, :), 'LineWidth', 1.2, 'Color', colors{ii});
end

for i=1:length(FP_Cell) 
rectangle('Position',[FP_Cell{i}(4)-5,FP_Cell{i}(16)-5,10,10],...
    'Curvature',[1,1],...
    'FaceColor','black'), hold on
daspect([1,1,1])
end  
    
set(gca, 'FontSize', 14)
xlabel('SNAI1', 'FontSize', 14); 
ylabel('PPAR', 'FontSize', 14);
pause(1); 
end