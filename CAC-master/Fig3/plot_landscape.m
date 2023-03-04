

%PLOT5ilandscape
a=ones(1,num);
%d1=12;d2=14;a(d1)=0.8;a(d2)=20;cn=1;
%d1=5;d2=16;a(d1)=5;a(d2)=1;cn=1;
d1=4;d2=16;a(d1)=15;a(d2)=5;cn=2;
%d1=4;d2=5;a(d1)=15;a(d2)=5;cn=3;
titl{4}={'xy=516','xy=416','xy=45'};
clen=1000;label={'ZEB1','Oct4','Mdm2','Snai1','miR145','miR200','miR34','p53','RKIP','Let7','Lin28','Bach1','MEK','ERK','CEBP','PPAR'};
weight=zeros(1,ssn);x_ss=zeros(num+num*num,ssn);cov_ss=zeros(ssn,num,num);sig_ss=zeros(num,ssn);
num=16;temp=1;P1=[];ds=max(YJ(1:num,:),[],2);dt=ds./clen.*3;
gx=clen;gy=clen;pt=zeros(gx,gy);
NN=size(YJ);
for ig=1:NN(2)
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

for j=1:16000
    for i=1:num
      P1(j,i,:)=exp(-power(j*dt(i)-x_ss(i,:),2)./2./sig_ss(i)/a(i))./sqrt(2.*pi.*sig_ss(i).*a(i));
    end
end
 
Pt=sum(P1,1);Pd=P1./Pt;
%G1
  for i=1:1:gx
  for j=1:1:gy
  for k=1:temp-1
  pt(i,j)=pt(i,j)+wei(1,k).*Pd(i,d1,k).*Pd(j,d2,k);
if pt(i,j)<=1e-75 %%%%solve probability zero
   pt(i,j)=1e-75;
end
   end
   end
  end
figure(1)
U15=-log(pt)';
[x1,x2] = meshgrid(0:dt(d1):(gx-1)*dt(d1), 0:dt(d2):(gy-1)*dt(d2));  
set(gcf,'outerposition', [100 80 750 600]);kx=1:1:gx;ky=1:1:gy;
surf(x1(ky,kx),x2(ky,kx),U15(ky,kx));shading interp;set(gca,'FontSize',16);
%xlabel('Snai1','FontSize',16),ylabel('miR34','FontSize',16),zlabel('U','FontSize',16),title('Fig.2','FontSize',16);
ylabel(label{d2},'FontSize',16),xlabel(label{d1},'FontSize',16),zlabel('U','FontSize',16),title('Fig.2','FontSize',16);
%axis([0 550 0 200 -45 200]);
%axis([0 (gx-1)*dt(d1) (gy-1)*dt(d2)*0.49 (gy-1)*dt(d2)*0.61 0 180]);

axis([0 (gx-1)*dt(d1)*0.7 0 (gy-1)*dt(d2)*0.7 0 180]);

%plotfig2,ad1=0.1,ad2=2;
%axis([0 8000 0 150 0 180]);
axis([0 1600 0 4600 0 180]);
%axis([0 160 0 150 0 180]);
view(281,87)  
hold on
%hash_ts=strcat('³ÌÐò\drugtreat\altercanshu\',hash_tr,titl{4}(l4),'.fig');
%hgsave(hash_ts{1}); 
%}
%
for i=1:ssn
    for j=1:ssn
        if i~=j&&~(i==2&&j~=1)&&~(i~=1&&j==2)
            hold on
           x_p=ycell{i,j}(d2,:);
y_p=ycell{i,j}(d1,:);
z_p=griddata(x1(ky,kx),x2(ky,kx),U15(kx,ky),y_p,x_p);
if i<j
plot3(y_p,x_p,z_p+10,'Color','m','LineWidth',2);
elseif i>j
plot3(y_p,x_p,z_p+10,'w','LineWidth',2);        
 end
%scatter3(y_p,x_p,z_p+10)
        end
    end
end
%}




