

weight=zeros(1,ssn);
x_ss=zeros(num+num*num,ssn);
cov_ss=zeros(ssn,num,num);
sig_ss=zeros(num,ssn);
%a=1.8;
%N=50;
temp=1;
for ig=1:ii
        if ~ismember(round(YJ(1,ig)),round(x_ss(1,:)))
         x_ss(:,temp)=YJ(:,ig);
       cov_ss(temp,:,:) = reshape(YJ(num+1:end,ig), num, num);
       sig_ss(:,temp)=diag(squeeze(cov_ss(temp,:,:)));
       temp=temp+1;
        end
        weight(1,find(x_ss(1,:)==YJ(1,ig)))= weight(1,find(x_ss(1,:)==YJ(1,ig)))+1;
end
wei=weight./sum(weight);
%{
for i=1:num
x_ss(i,:)=x_ss(i,:)./max(x_ss(i,:));
end
%}
  SS=x_ss(1:num,:);


params={}; 
params.N=80;%59 %The number of points in the minimum action path.6~15 <20
params.TMax = 1;%Ö®Ç°ÊÇ1,10  %The time range of the minimum action path. %100
%Larger values can be more accurate, but can also lead to instabilities. 
params.Dimension = size(SS, 1); %The dimension of the system. 
params.c = 1e12; %The remeshing parameter, c. Larger values give greater remeshing. 
params.MaxIter = 300; %The number of optimization steps to run between remeshing steps. 
params.K= 2; %The number of total remeshing steps to do; 
params.Remesh = 1; %Turn the remesh on or off. If off, the algorithm will be default perform K*MaxIter itereations. 
params.q = 3; %The q parameter from the original paconstraints can be set of the path as well. This is not done in our per. Is a measure of path smoothness. 
params.ub = []; %If desired, manuscript. 
params.lb = []; %As above, lower bounds can be set on the path as well. This is also not done at all in our manuscript. 
params.PhiInit =[]; %The default intial path between two states is a straight line. This can be modified here. 
%params.ObjGrad = Parameters.JacobianPresent;%We require a Jacobian in this implementation. 
params.ObjGrad=0;

Func=@(x)Force5i(1,x,y);
%Func=@(x)Force5jj(1,x,y)
%dFunc=@(x)jacobi(1,x,y);
dFunc=[];

feasi_n=size(SS,2);


action=zeros(feasi_n,feasi_n);
ycell=cell(feasi_n,feasi_n); %%path
for i=1:feasi_n
    for j=1:feasi_n
        if i==j
            action(i,j) = Inf;
        else
            Spot1 = i;
            Spot2 = j;%initial and end point
            params.init = SS(:, [Spot1,Spot2]);
            [ActionVal, Path] = AMAM(params, Func, dFunc);
            action(i,j) = ActionVal;
            ycell(i,j) = {[SS(:,i) Path SS(:,j)]};
        end
        (i-1)*feasi_n+j
    end
end
action

