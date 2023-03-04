
num=16;

YY=[];YJ=[];N=40;nbatch=1;T=4000;h=0.005;a=10;n=2;d=0.3;nstep=T/h;l=800;m=10;u=0.2;
load('stable.mat');
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
%for n=1:20
out=0;

for ii=1:N
x=randi(3000,num+num*num,nbatch);
ii
if ii<6
x(1:num)=stable(1:num,ii);
end
ff=zeros(num+num*num, nbatch);
%
for j=1:nstep
          ff=cov_i(nstep, x, y);
          res=h*ff;
          if norm(res(1:16))<1e-7
              j
              break
          end
          x=x+res;
          
end
%}
 YJ(:,ii,:)=x;
 YJ1=reshape(YJ,num+num*num,[]);
 for g=1:num
 ssnt(g)=length(unique(round(YJ1(g,:))));    
 end
 ssn=ssnt(1)
 if ssn>4
     break
 end
ii
end

YJ=reshape(YJ,num+num*num,[]);

ssnt=[];
for g=1:num
 ssnt(g)=length(unique(round(YJ(g,:))))
end
ssn=ssnt(1);
cal_path5i;
plot_landscape;
