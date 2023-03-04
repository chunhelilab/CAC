Var_label={
'ZEB1/2 (x_{1})', 'OCT4 (x_{2})', 'MDM2 (x_{3})', 'SNAIL (x_{4})', 'miR-145 (x_{5})', 'miR-200 (x_{6})',...
'miR-34 (x_{7})', 'P53 (x_{8})', 'RKIP (x_{9})', 'Let7 (x_{10})', 'LIN28 (x_{11})',...
'BACH1 (x_{12})', 'MEK1/2 (x_{13})', 'ERK1/2 (x_{14})', 'CEBP_{/alpha} (x_{15})', 'PPAR_{/gamma} (x_{16})'};

num=16;TT={};
YY={};YJ=[];N=1000;nbatch=1;h=0.005;T=2000;nstep=T/h;
a=10;n=2;n2=4;d=0.3;l=800;m=10;u=0.2;m2=15; 
y=[ a  a  a  a	d  d  d  d	u	24,...      %1-10
    m  m  m  m	m  m  m  4	n	4,...     %11-20
    n  n  n  n	n  n  l  l	1000 	l,...  %21-30
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
    n n n n n n n n n n,...                             %141-150
    l,l,l,l,l,l,l,l,l,l,...         %151-160
    l,l,l,l,l,l,l,l,l,l,...         %161-170
    l,l,l,l,20,l,l,l,l,l,...         %171-180
    l,l,l,l,u,n,l,m,n];                                       %181-190
load stable.mat
st=1/h;  
for ii=1:N
if ii<=4
    x=stable(:,ii);
else
x=randi(3000,num,nbatch);
end
ff=zeros(num,nbatch);
for j=1:nstep
 ff=Force5i(nstep,x,y);
          res=h.*ff;
          if norm(res(1:16),1)<1e-7
              j
              break
          end
          x=x+res;
  if mod(j,st)==0
    YY{ii}(:,j/st,:)=x(:,:);
    TT{ii}(j/st)=j*h;
 end
end
 YJ(:,ii,:)=x;
  %if mod(ii,10)==0 
  YJ1=reshape(round(YJ),num,[]);
  for g=1:num
 ssnt(g)=length(unique(round(YJ1(g,:))));    
 end
 ssn=ssnt(1) 
 
 if ssn>3
     break
 end
  %end
 ii
end

res=sort(unique(round(YJ(1,:))),'ascend');
idx1=find(YJ1(1,:)==res(1));
idx2=find(YJ1(1,:)==res(2));
idx3=find(YJ1(1,:)==res(3));
idx4=find(YJ1(1,:)==res(4));
stable=[YJ(:,idx1(1)),YJ(:,idx2(1)),YJ(:,idx4(1)),YJ(:,idx3(1));];
save('y_4ss','y')
save('stable_4ss','stable')
Optimization
OLAC_param_plot