



Var_label={
'ZEB1/2 (x_{1})', 'OCT4 (x_{2})', 'MDM2 (x_{3})', 'SNAIL (x_{4})', 'miR-145 (x_{5})', 'miR-200 (x_{6})',...
'miR-34 (x_{7})', 'P53 (x_{8})', 'RKIP (x_{9})', 'Let7 (x_{10})', 'LIN28 (x_{11})',...
'BACH1 (x_{12})', 'MEK1/2 (x_{13})', 'ERK1/2 (x_{14})', 'CEBP_{/alpha} (x_{15})', 'PPAR_{/gamma} (x_{16})'};
clc;
clear;

result=[];
while size(result,1)==0
    result=search_param(result);
    if size(result,1)>0
        break
    end
end



function result=search_param(result)
param_list={'a','d','m','u','n','L','TGF','MEKi','Rosi'};
range_partitions=[0 20;
                  0 2;
                  0 20;
                  0 1;
                  2 4;
                  0 2000;
                  1 6;
                  0.01 1.2;
                  0 30];
partitions_length=[1,0.1,1,0.1,1,100,0.1,0.01,0.1];
%a=10;n=2;d=0.3;L=800;m=10;u=0.2;
param=initial_conditions(param_list,range_partitions,partitions_length);
y=set_y(param);
done=0;stepped=5;
ssn=runode(y);

 for index=1:6
        index
        done=0;stp=0;over=0;notgoback=0;
    while stp<=stepped
        new_param=param;
       if ~done
        new_param(index)=param(index)+partitions_length(index);
       else
        new_param(index)=param(index)-partitions_length(index);
       end
       if new_param(index)>=range_partitions(index,2) || new_param(index)<=range_partitions(index,1)
           stp=stp+1;
            continue
       end           
        newy=set_y(new_param);
        new_ssn=runode(newy);
        if new_ssn>2
            result=[result;new_param];
            break
        end
        if new_ssn>ssn
            notgoback=1;
            stp=0;
            param=new_param;
            ssn=new_ssn;
            if ~done
            done=~done;
            end
        else
            stp=stp+1;
        end
      if stp==stepped && notgoback
          break
      end
    end
end
end
function y=set_y(p)
a=p(1);d=p(2);m=p(3);u=p(4);n=p(5);L=p(6);
TGF=p(7);MEKi=p(8);Rosi=p(9);
y=[ a  a  a  a	d  d  d  d	 u	 m,...      %1-10
    m  m  m  m	m  m  m  n	 n	 n,...     %11-20
    n  n  n  n	n  n  L  L	L L,...  %21-30
    L L L  L L  TGF L MEKi  Rosi  m,... %31-40
    L n a L n L L L L L,...         %41-50
    a a a a a a a a a a,...         %51-60
    a a d d d d d d d d,...         %61-70
    d d d d m m m m u u,...         %71-80
    m m m m u m m m m m,...         %81-90
    m m m m m m m m u m,...         %91-100
    m m u u m m u m m m,...         %101-110                 
    m m n n n n n n n n,...         %111-120
    n n n n n n n n n n,...         %121-130
    n n n n n n n n n n,...         %131-140
    n n n n n n n n n n,...         %141-150
    L,L,L,L,L,L,L,L,L,L,...         %151-160
    L,L,L,L,L,L,L,L,L,L,...         %161-170
    L,L,L,L,a,L,L,L,L,L,...         %171-180
    L,L,L,L,u,n,L,m,n];             %181-190
end

function param=initial_conditions(param_list,range_partitions,partitions_length)
res=range_partitions(:,1)+round(unifrnd(range_partitions(:,1),range_partitions(:,2))./partitions_length').*partitions_length';
param=res';
end

function ssn=runode(y)
h=0.005;num=16;TT={};YY={};YJ=[];N=300;nbatch=1;T=2000;nstep=T/h;st=1/h; 
for ii=1:N
x=randi(3000,num,nbatch);
%YY(:,ii,:)=round(x);
ff=zeros(num,nbatch);
for j=1:nstep
 ff=BFS_model(nstep,x,y);
          res=h.*ff;
          if norm(res(1:16),1)<1e-7
              %j;
              break
          end
          x=x+res;
  if mod(j,st)==0
    YY{ii}(:,j/st,:)=x(:,:);
    TT{ii}(j/st)=j*h;
 end
end
 YJ(:,ii,:)=round(x);
  %if mod(ii,10)==0 
  YJ1=reshape(YJ,num,[]);
  for g=1:num
 ssnt(g)=length(unique(round(YJ1(g,:))));    
 end
 ssn=max(ssnt);
  if sum(sum(1*(isnan(YJ))))>0
      ssn=0
      break
  end
 
 if ssn>2
     YJ
     break
 end
 if mod(ii,150)==0
 ii
 ssn
 end
end
end