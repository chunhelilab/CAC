
Param_Label={'g_{MEK}','g_{ERK}','g_{CEBP}','g_{PPAR}','k_{MEK}','k_{ERK}','k_{CEBP}','k_{PPAR}','lambda_{EM}','lambda_{ME}','lambda_{EE}','lambda_{MC}','lambda_{PC}','lambda_{EC}','lambda_{MP}','lambda_{CP}','lambda_{EP}','n_{EM}','n_{ME}','n_{EE}','n_{MC}','n_{PC}','n_{EC}','n_{MP}','n_{CP}','n_{EP}','MEK_{MC}','MEK_{MP}','MEK0_{ME}','ERK0_{EC}','ERK0_{EP}','CEBP0_{CP}','PPAR_{PC}','ERK0_{EM}','ERK0_{EE}','TGF_\beta','rkip0_{RM}','MEK_i','Ros_i','lambda_{EEi}','ERK0_{EEi}','n_{EEi}','lambda_{MM}','MEK0_{MM}','n_{MM}','MEK0_{Ml}','rkip0_{RS}','oct40_{ol}','bach0_{br}','bach0_{bb}','g_{Z}','g_{o}','g_{M}','g_{S}','g_{m1}','g_{m2}','g_{m3}','g_{p}','g_{rk}','g_{let}','g_{LIN}','g_{b}','k_{Z}','k_{o}','k_{M}','k_{S}','k_{m1}','k_{m2}','k_{m3}','k_{p53}','k_{rk}','k_{let}','k_{LIN}','k_{b}','lambda_{sZ}','lambda_{ZZ}','lambda_{m2Z}','lambda_{m1Z}','lambda_{pO}','lambda_{oo}','lambda_{lo}','lambda_{m1o}','lambda_{pM}','lambda_{m1M}','lambda_{SS}','lambda_{m3s}','lambda_{ls}','lambda_{om1}','lambda_{Zm1}','lambda_{pm1}','lambda_{om2}','lambda_{pm2}','lambda_{Zm2}','lambda_{sm2}','lambda_{sm3}','lambda_{Zm3}','lambda_{pm3}','lambda_{Mp}','lambda_{br}','lambda_{sr}','lambda_{RL}','lambda_{llet}','lambda_{let}','lambda_{LIN}','lambda_{llin}','lambda_{m2le}','lambda_{bb}','lambda_{lb}','lambda_{RM}','lambda_{Ml}','lambda_{RS}','lambda_{ol}','n_{sZ}','n_{ZZ}','n_{m2Z}','n_{m1Z}','n_{pO}','n_{oo}','n_{lo}','n_{m1o}','n_{pM}','n_{m1M}','n_{SS}','n_{m3s}','n_{ls}','n_{om1}','n_{Zm1}','n_{pm1}','n_{om2}','n_{pm2}','n_{Zm2}','n_{sm2}','n_{sm3}','n_{Zm3}','n_{pm3}','n_{Mp}','n_{br}','n_{sr}','n_{RL}','n_{llet}','n_{let}','n_{LIN}','n_{llin}','n_{m2le}','n_{bb}','n_{lb}','n_{RM}','n_{Ml}','n_{RS}','n_{ol}','Z0_{ZZ}','Z0_{Zm1}','Z0_{Zm2}','Z0_{Zm3}','oct40_{oo}','oct40_{om1}','oct40_{om2}','MD0_{Mp}','snail0_{sZ}','snail0_{SS}','snail0_{sm2}','snail0_{sm3}','snail0_{sr}','miR1450_{m1Z}','miR1450_{m1o}','miR1450_{m1M}','miR2000_{m2Z}','miR2000_{m2le}','miR340_{m3s}','p0_{pO}','p0_{pM}','p0_{pm1}','p0_{pm2}','p0_{pm3}','P530}','rkip0_{RL}','let0_{ls}','let0_{let}','let0_{ll}','let0_{lb}','lin0_{lo}','lin0_{llet}','lin0_{LIN}','rkip0_{rli}','lambda_{rli}','n_{rli}','miR2000_{m2li}','lambda_{m2li}','n_{m2li}'};
PL_use1={};adim=189;
a=10;n=2;d=0.3;l=800;m=10;u=0.2;n2=4;m2=15;

load y_4ss.mat
load OutPut_4ss.mat
id=1
%for id=1:7
%y_path=strcat('good_param/y_',num2str(id),'.mat');
%OutPut_path=strcat('good_param/OutPut_',num2str(id),'.mat');
%save_path=strcat('good_param_result/result_',num2str(id),'.png');
%load(y_path)
%load(OutPut_path)

result_figure=figure(100-id)
nn=189;
Param_Label={'g_{MEK}','g_{ERK}','g_{CEBP}','g_{PPAR}','k_{MEK}','k_{ERK}','k_{CEBP}','k_{PPAR}','\lambda_{EM}','\lambda_{ME}','\lambda_{EE}','\lambda_{MC}','\lambda_{PC}','\lambda_{EC}','\lambda_{MP}','\lambda_{CP}','\lambda_{EP}','n_{EM}','n_{ME}','n_{EE}','n_{MC}','n_{PC}','n_{EC}','n_{MP}','n_{CP}','n_{EP}','MEK_{MC}','MEK_{MP}','MEK0_{ME}','ERK0_{EC}','ERK0_{EP}','CEBP0_{CP}','PPAR_{PC}','ERK0_{EM}','ERK0_{EE}','TGF{\beta}','rkip0_{RM}','MEKi','Rosi','\lambda_{EEi}','ERK0_{EEi}','n_{EEi}','\lambda_{MM}','MEK0_{MM}','n_{MM}','MEK0_{Ml}','rkip0_{RS}','oct40_{ol}','bach0_{BR}','bach0_{BB}','g_{Z}','g_{o}','g_{M}','g_{S}','g_{m1}','g_{m2}','g_{m3}','g_{p}','g_{RKIP}','g_{let}','g_{LIN}','g_{BACH1}','k_{ZEB}','k_{OCT4}','k_{MDM2}','k_{S}','k_{m1}','k_{m2}','k_{m3}','k_{p}','k_{RKIP}','k_{let}','k_{LIN}','k_{BACH1}','\lambda_{sZ}','\lambda_{ZZ}','\lambda_{m2Z}','\lambda_{m1Z}','\lambda_{pO}','\lambda_{oo}','\lambda_{lo}','\lambda_{m1o}','\lambda_{pM}','\lambda_{m1M}','\lambda_{SS}','\lambda_{m3s}','\lambda_{ls}','\lambda_{om1}','\lambda_{Zm1}','\lambda_{pm1}','\lambda_{om2}','\lambda_{pm2}','\lambda_{Zm2}','\lambda_{sm2}','\lambda_{sm3}','\lambda_{Zm3}','\lambda_{pm3}','\lambda_{Mp}','\lambda_{BR}','\lambda_{sr}','\lambda_{RL}','\lambda_{llet}','\lambda_{let}','\lambda_{LIN}','\lambda_{llin}','\lambda_{m2le}','\lambda_{BB}','\lambda_{LB}','\lambda_{RM}','\lambda_{Ml}','\lambda_{RS}','\lambda_{ol}','n_{sZ}','n_{ZZ}','n_{m2Z}','n_{m1Z}','n_{pO}','n_{oo}','n_{lo}','n_{m1o}','n_{pM}','n_{m1M}','n_{SS}','n_{m3s}','n_{ls}','n_{om1}','n_{Zm1}','n_{pm1}','n_{om2}','n_{pm2}','n_{Zm2}','n_{sm2}','n_{sm3}','n_{Zm3}','n_{pm3}','n_{Mp}','n_{BR}','n_{sr}','n_{RL}','n_{llet}','n_{let}','n_{LIN}','n_{llin}','n_{m2le}','n_{BB}','n_{LB}','n_{RM}','n_{Ml}','n_{RS}','n_{ol}','Z0_{ZZ}','Z0_{Zm1}','Z0_{Zm2}','Z0_{Zm3}','oct40_{oo}','oct40_{om1}','oct40_{om2}','MD0_{Mp}','snail0_{sZ}','snail0_{SS}','snail0_{sm2}','snail0_{sm3}','snail0_{sr}','miR1450_{m1Z}','miR1450_{m1o}','miR1450_{m1M}','miR2000_{m2Z}','miR2000_{m2le}','miR340_{m3s}','p0_{pO}','p0_{pM}','p0_{pm1}','p0_{pm2}','p0_{pm3}','P530','rkip0_{RL}','let0_{ls}','let0_{let}','let0_{ll}','let0_{LB}','lin0_{lo}','lin0_{llet}','lin0_{LIN}','rkip0_{rli}','\lambda_{rli}','n_{rli}','miR2000_{m2li}','\lambda_{m2li}','n_{m2li}'};
PL_use={};
y0=y;
param=OutPut.Params(:);
[ath,indexath]=sort((param'-y0)./y0,'descend');
%[ath,indexath1]=sort(abs(param'-y0)./y0,'descend');
for i=1:189
    if (y0(indexath(i))-param(indexath(i)))>0
        ath(i)=ath(i)*(1);
    else
        ath(i)=ath(i)*(1);
    end
end

ath1=[ath(1:5)';ath(184:189)'];

for i=1:5
PL_use(i)=Param_Label(indexath(i));
end
%
for i=6:10
PL_use(i)=Param_Label(indexath(i-5+184));
end
%}
bh=bar(ath1,0.8)
xlim([0,11-0.5])
bh(1).FaceColor=[1 0.3 0.1]
hold on 


x1=set(gca, 'XTick', 1:1:10, 'XTickLabel', PL_use,'FontSize',12,'FontWeight','bold','FontName','Arial');
set(gca,'XTickLabelRotation',90)


set(gca,'Box','off')
caxis([0 1])
%xlabel('FontSize',14,'FontWeight','bold','FontName','Arial');
% ylabel('FontSize',14,'FontWeight','bold','FontName','Arial');

[M,N]=size(ath1);
set(gca,'Clim',[0 1])
set(result_figure, 'unit', 'normalized', 'position', [0.1,0.05,0.35,0.3])

C=zeros(2,3);
C(1,1)=1;
C(1,2)=0;
C(1,3)=1;
C(2,1)=0;
C(2,2)=1;
C(2,3)=1;
colormap(C);
%saveas(100-id, save_path); 
%end