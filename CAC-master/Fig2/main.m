%This is an implementation of the DRL(dimension reduction of landscape)
%method on CAC model.
cycle_index=800;  
%% The number of random initial conditions to the ODEs to be solved
global num;num=16;YY=[];YJ=[];N=40;nbatch=1;T=4000;
h=0.005;a=10;n=2;d=0.3;nstep=T/h;l=800;m=10;u=0.2;
y=[ a  a  a  a	d  d  d  d	u	 24,...      %1-10
    m  m  m  m	m  m  m  4	n	 4,...     %11-20
    n  n  n  n	n  n  l  l	1000 l,...  %21-30
    l l l  l l  4.3 l 0.02  1.5   m,... %31-40
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
    n n n n n n n n n n,...                             %141-150
    l,l,l,l,l,l,l,l,l,l,...         %151-160
    l,l,l,l,l,l,l,l,l,l,...         %161-170
    l,l,l,l,20,l,l,l,l,l,...         %171-180
    l,l,l,l,u,n,l,m,n];                                       %181-190
%% The dimension of the system
global Num;
Num=16;
tic() %% Time
diffusion=30;  %% The diffusion coefficient 
%% Solve the ODEs, calculate the paths and actions;
[xx,sigma,n,ycell,action]=Solver(cycle_index,y,diffusion);
%xx=s_s;
index=size(n,1);  %% The number of the stable states
alpha=zeros(index,1);  %% The weight of the stable states
sigma0=cell(index,1);  %% The covariance of the Gaussian density function
mu=zeros(index,Num);  %% The mean value of the Gaussian density function

for i=1:index
   %The mean value of each stable state
   mu(i,:)=xx(n(i,1),:); 
   %The covariance of each stable state
   sigma0{i}=reshape(sigma(n(i,1),:),Num,Num)';  
   %The weight of each stable state
   alpha(i)=n(i,2)/sum(n(:,2)); 
end

%% DRL
%Calculate the mean value 
Mu=0;
for i=1:index
    Mu=Mu+alpha(i)*mu(i,:);
end
%Calculate the covariance
Sigma=-Mu'*Mu; 
for i=1:index
    Sigma=Sigma+alpha(i)*(sigma0{i}+mu(i,:)'*mu(i,:));
end

%Calculate the eigenvalues and eigenvectors of the covariance
[V,D] = eigs(Sigma,2);

if sign(V(:,1)'*[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]')<0
    V(:,1)=-V(:,1);
end
if sign(V(:,2)'*[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]')<0
    V(:,2)=-V(:,2);
end

%%%Calculate the covariance and mean value after dimension reduction
sigma0_pca=cell(index,1);
mu_pca=zeros(index,2);
for i=1:index
   mu_pca(i,:)=V'*mu(i,:)';
   sigma0_pca{i}=V'*sigma0{i}*V;
end
%
at=1
sigma0_pca{2}=sigma0_pca{1}*at
sigma0_pca{3}=sigma0_pca{1}*at
sigma0_pca{4}=sigma0_pca{1}*at
sigma0_pca{5}=sigma0_pca{1}*at
sigma0_pca{1}=sigma0_pca{1}*at

%}
%{
at=2;
sigma0_pca{1}=sigma0_pca{1}*at
sigma0_pca{2}=sigma0_pca{2}*at
sigma0_pca{3}=sigma0_pca{3}*at
sigma0_pca{4}=sigma0_pca{4}*at
sigma0_pca{5}=sigma0_pca{5}*at
%}
%{
sigma0_pca{1}=sigma0_pca{1}*1.3
sigma0_pca{2}=sigma0_pca{1}
sigma0_pca{3}=sigma0_pca{1}
sigma0_pca{4}=sigma0_pca{1}
sigma0_pca{5}=sigma0_pca{1}
%}
% plot the landscape
p1_max=15500;
p1_min=-4000;
p2_max=3500;
p2_min=-7000;
y_max=[p1_max,p2_max]; %% Range of the landscape
y_min=[p1_min,p2_min];
step=(y_max-y_min)/1250; %% Length of the step
[a1,a2]=meshgrid(y_min(1):step(1):y_max(1),y_min(2):step(2):y_max(2)); %% Grid
[s1,s2]=size(a1);
P=zeros(s1,s2);
z=zeros(s1,s2);
for kk=1:index
    sig=sigma0_pca{kk};
    x_wen=mu_pca(kk,:);
    for i=1:s1
        for j=1:s2
            z(i,j)=multivariate_normal_distribution([a1(i,j);a2(i,j)],x_wen',sig,2);  %% Normal distribution
        end
    end

    P=P+z*alpha(kk);
end
P=P/sum(sum(P));
surf(a1,a2,-log(max(P,10^-100)));   %% Plot landscape
shading interp
xlabel('PC1','FontSize',12)
ylabel('PC2','FontSize',12)
zlabel('U','FontSize',12)
ax=gca;
ax.ZGrid = 'on';
axis([p1_min p1_max p2_min p2_max 0 250])
set(gca,'FontSize',12);
for i=1:size(n,1)
    A(i)=floor((mu_pca(i,1)-y_min(1))/step(1))+1;
    B(i)=floor((mu_pca(i,2)-y_min(2))/step(2))+1;
end
 hold on
%{
%Plot the grid
for i=1:floor(size(a1,1)/4)
    plot3(a1(4*i-1,:),a2(4*i-1,:),-log(max(P(4*i-1,:),10^-100)),'Color',[0.4 0.4 0.4],'LineWidth',0.01);
end
for i=1:floor(size(a1,2)/4)
    plot3(a1(:,4*i-1),a2(:,4*i-1),-log(max(P(:,4*i-1),10^-100)),'Color',[0.4 0.4 0.4],'LineWidth',0.01);
end
%}
%
%Calculate the paths after dimension reduction
k=size(ycell);


y12=V'*ycell{5,4}
hold on
z3path=griddata(a1,a2,-log(max(P,10^-100)),y12(1,:),y12(2,:));
plot3(y12(1,:),y12(2,:),z3path+10,'Color',[0.85,0.43,0.83],'LineWidth',2);
view([-25,75])   

y12=V'*ycell{4,2}
hold on
z3path=griddata(a1,a2,-log(max(P,10^-100)),y12(1,:),y12(2,:));
plot3(y12(1,:),y12(2,:),z3path+10,'Color',[0.85,0.43,0.83],'LineWidth',2);
view([-25,75])   

y12=V'*ycell{2,1}
hold on
z3path=griddata(a1,a2,-log(max(P,10^-100)),y12(1,:),y12(2,:));
plot3(y12(1,:),y12(2,:),z3path+10,'Color',[0.85,0.43,0.83],'LineWidth',2);
view([-25,75])   



y12=V'*ycell{1,3}
hold on
z3path=griddata(a1,a2,-log(max(P,10^-100)),y12(1,:),y12(2,:));
plot3(y12(1,:),y12(2,:),z3path+10,'y','LineWidth',2);
view([-25,75])   

y12=V'*ycell{5,3}
hold on
z3path=griddata(a1,a2,-log(max(P,10^-100)),y12(1,:),y12(2,:));
plot3(y12(1,:),y12(2,:),z3path+10,'y','LineWidth',2);
view([-25,75])   





%{
for i=1:k(2)
    for j=1:k(2)
        if i-j==2 && (i==1||i==3||i==5)
  y12=V'*ycell{i,j}
hold on
z3path=griddata(a1,a2,-log(max(P,10^-100)),y12(1,:),y12(2,:));
plot3(y12(1,:),y12(2,:),z3path+10,'r','LineWidth',2);
view([-25,75])      
        end
         if j-i==2 && (i==1||i==3||i==5)
  y12=V'*ycell{i,j}
hold on
z3path=griddata(a1,a2,-log(max(P,10^-100)),y12(1,:),y12(2,:));
plot3(y12(1,:),y12(2,:),z3path+10,'g','LineWidth',2);
view([-25,75])      
        end
               
        
        
        if i-j==1
y12=V'*ycell{i,j}
hold on
z3path=griddata(a1,a2,-log(max(P,10^-100)),y12(1,:),y12(2,:));
plot3(y12(1,:),y12(2,:),z3path+10,'w','LineWidth',2);
view([-25,75])
        elseif j-i==1
y21=V'*ycell{i,j}
hold on 
z3path=griddata(a1,a2,-log(max(P,10^-100)),y21(1,:),y21(2,:));
plot3(y21(1,:),y21(2,:),z3path+10,'Color',[0.85,0.43,0.83],'LineWidth',2);
view([-25 75])
set(gcf,'outerposition', [100 100 800 650]);
toc()
        end
    end
end
%}