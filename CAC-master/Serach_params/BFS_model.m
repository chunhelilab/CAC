
function df=BFS_model(t,x,y)
%%Components of the model
Z=x(1,:);oct4=x(2,:);MD=x(3,:);snail=x(4,:);miR145=x(5,:);miR200=x(6,:);miR34=x(7,:);p53=x(8,:);
rkip=x(9,:);let=x(10,:);lin=x(11,:);bach=x(12,:);MEK=x(13,:);
ERK=x(14,:);CEBP=x(15,:);PPAR=x(16,:);

%%Parameters of the model
g_MEK=y(1);g_ERK=y(2);g_CEBP=y(3);g_PPAR=y(4);k_MEK=y(5);k_ERK=y(6);k_CEBP=y(7);k_PPAR=y(8);lambda_EM=y(9);lambda_ME=y(10);
lambda_EE=y(11);lambda_MC=y(12);lambda_PC=y(13);lambda_EC=y(14);lambda_MP=y(15);lambda_CP=y(16);lambda_EP=y(17);n_EM=y(18);n_ME=y(19);n_EE=y(20);
n_MC=y(21);n_PC=y(22);n_EC=y(23); n_MP=y(24);n_CP=y(25);n_EP=y(26);MEK_MC=y(27);MEK_MP=y(28);MEK0_ME=y(29);ERK0_EC=y(30);
ERK0_EP=y(31);CEBP0_CP=y(32);PPAR_PC=y(33);ERK0_EM=y(34);ERK0_EE=y(35);TGF=y(36);rkip0_rM=y(37);MEKi=y(38);Rosi=y(39);lambda_EEi=y(40);
ERK0_EEi=y(41);n_EEi=y(42);lambda_MM=y(43);MEK0_MM=y(44);n_MM=y(45);MEK0_Ml=y(46);rkip0_rs=y(47);oct40_ol=y(48);bach0_br=y(49);bach0_bb=y(50);

g_Z=y(51);g_o=y(52);g_M=y(53);g_s=y(54);g_m1=y(55);g_m2=y(56);g_m3=y(57);g_p=y(58);g_rk=y(59);g_let=y(60);
g_lin=y(61);g_b=y(62);k_Z=y(63);k_o=y(64);k_M=y(65);k_s=y(66);k_m1=y(67);k_m2=y(68);k_m3=y(69);k_p=y(70);
k_rk=y(71);k_let=y(72);k_lin=y(73);k_b=y(74);lambda_sZ=y(75);lambda_ZZ=y(76);lambda_m2Z=y(77);lambda_m1Z=y(78);lambda_po=y(79);lambda_oo=y(80);
lambda_lo=y(81);lambda_m1o=y(82);lambda_pM=y(83);lambda_m1M=y(84);lambda_ss=y(85);lambda_m3s=y(86);lambda_ls=y(87);lambda_om1=y(88);lambda_Zm1=y(89);lambda_pm1=y(90);
lambda_om2=y(91);lambda_pm2=y(92);lambda_Zm2=y(93);lambda_sm2=y(94);lambda_sm3=y(95);lambda_Zm3=y(96);lambda_pm3=y(97);lambda_Mp=y(98);lambda_br=y(99);lambda_sr=y(100);
lambda_rl=y(101);lambda_llet=y(102);lambda_let=y(103);lambda_lin=y(104);lambda_llin=y(105);lambda_m2le=y(106);lambda_bb=y(107);lambda_lb=y(108);lambda_rM=y(109);lambda_Ml=y(110);
lambda_rs=y(111);lambda_ol=y(112);n_sZ=y(113);n_ZZ=y(114);n_m2Z=y(115);n_m1Z=y(116);n_po=y(117);n_oo=y(118);n_lo=y(119);n_m1o=y(120);
n_pM=y(121);n_m1M=y(122);n_ss=y(123);n_m3s=y(124);n_ls=y(125);n_om1=y(126);n_Zm1=y(127);n_pm1=y(128);n_om2=y(129);n_pm2=y(130);
n_Zm2=y(131);n_sm2=y(132);n_sm3=y(133);n_Zm3=y(134);n_pm3=y(135);n_Mp=y(136);n_br=y(137);n_sr=y(138);n_rl=y(139);n_llet=y(140);
n_let=y(141);n_lin=y(142);n_llin=y(143);n_m2le=y(144);n_bb=y(145);n_lb=y(146);n_rM=y(147);n_Ml=y(148);n_rs=y(149);n_ol=y(150);
Z0_ZZ=y(151);Z0_Zm1=y(152);Z0_Zm2=y(153);Z0_Zm3=y(154);oct40_oo=y(155);oct40_om1=y(156);oct40_om2=y(157);MD0_Mp=y(158);snail0_sZ=y(159);snail0_ss=y(160);
snail0_sm2=y(161);snail0_sm3=y(162);snail0_sr=y(163);miR1450_m1Z=y(164);miR1450_m1o=y(165);miR1450_m1M=y(166);miR2000_m2Z=y(167);miR2000_m2le=y(168);miR340_m3s=y(169);p0_po=y(170);
p0_pM=y(171);p0_pm1=y(172);p0_pm2=y(173);p0_pm3=y(174);P530=y(175);rkip0_rl=y(176);let0_ls=y(177);let0_let=y(178);let0_ll =y(179);let0_lb=y(180);
lin0_lo=y(181);lin0_llet=y(182);lin0_lin=y(183);rkip0_rli=y(184);lambda_rli=y(185);n_rli=y(186);miR2000_m2li=y(187);lambda_m2li=y(188);n_m2li=y(189);

df(1,:)=g_Z.*Hs(snail,snail0_sZ,lambda_sZ,n_sZ).*TGF.*Hs(Z,Z0_ZZ,lambda_ZZ,n_ZZ)-k_Z.*Z.*Hs(miR200,miR2000_m2Z,lambda_m2Z,n_m2Z).*Hs(miR145,miR1450_m1Z,lambda_m1Z,n_m1Z);
df(2,:)=g_o.*Hs(p53,p0_po,lambda_po,n_po).*Hs(lin,lin0_lo,lambda_lo,n_lo)-k_o.*oct4.*Hs(miR145,miR1450_m1o,lambda_m1o,n_m1o).*Hs(oct4,oct40_oo,lambda_oo,n_oo);
df(3,:)=g_M.*Hs(p53,p0_pM,lambda_pM,n_pM)-k_M.*MD.*Hs(miR145,miR1450_m1M,lambda_m1M,n_m1M);
df(4,:)=g_s.*Hs(snail,snail0_ss,lambda_ss,n_ss).*TGF-k_s.*snail.*Hs(miR34,miR340_m3s,lambda_m3s,n_m3s).*Hs(let,let0_ls,lambda_ls,n_ls).*Hs(rkip,rkip0_rs,lambda_rs,n_rs);
df(5,:)=g_m1.*Hs(p53,p0_pm1,lambda_pm1,n_pm1)-k_m1.*miR145.*Hs(Z,Z0_Zm1,lambda_Zm1,n_Zm1).*Hs(oct4,oct40_om1,lambda_om1,n_om1);
df(6,:)=g_m2.*Hs(p53,p0_pm2,lambda_pm2,n_pm2).*Hs(oct4,oct40_om2,lambda_om2,n_om2)-k_m2.*miR200.*Hs(Z,Z0_Zm2,lambda_Zm2,n_Zm2).*Hs(snail,snail0_sm2,lambda_sm2,n_sm2);
df(7,:)=g_m3.*Hs(p53,p0_pm3,lambda_pm3,n_pm3)-k_m3.*miR34.*Hs(Z,Z0_Zm3,lambda_Zm3,n_Zm3).*Hs(snail,snail0_sm3,lambda_sm3,n_sm3);
df(8,:)=g_p.*P530-k_p.*p53.*Hs(MD,MD0_Mp,lambda_Mp,n_Mp);
df(9,:)=g_rk.*Hs(bach,bach0_br,lambda_br,n_br)-k_rk.*rkip.*Hs(snail,snail0_sr,lambda_sr,n_sr);
df(10,:)=g_let.*Hs(miR200,miR2000_m2le,lambda_m2le,n_m2le).*Hs(rkip,rkip0_rl,lambda_rl,n_rl)-k_let.*let.*Hs(lin,lin0_llet,lambda_llet,n_llet).*Hs(let,let0_let,lambda_let,n_let);
df(11,:)=g_lin.*Hs(MEK,MEK0_Ml,lambda_Ml,n_Ml).*Hs(oct4,oct40_ol,lambda_ol,n_ol).*Hs(rkip,rkip0_rli,lambda_rli,n_rli)-k_lin.*lin.*Hs(let,let0_ll,lambda_llin,n_llin).*Hs(lin,lin0_lin,lambda_lin,n_lin).*Hs(miR200,miR2000_m2li,lambda_m2li,n_m2li);
df(12,:)=g_b.*Hs(bach,bach0_bb,lambda_bb,n_bb)-k_b.*bach.*Hs(let,let0_lb,lambda_lb,n_lb);
df(13,:)=g_MEK.*TGF.*Hs(MEK,MEK0_MM,lambda_MM,n_MM)-k_MEK.*MEK.*Hs(rkip,rkip0_rM,lambda_rM,n_rM).*MEKi.*Hs(ERK,ERK0_EM,lambda_EM,n_EM);
df(14,:)=g_ERK.*Hs(MEK,MEK0_ME,lambda_ME,n_ME).*Hs(ERK,ERK0_EE,lambda_EE,n_EE)-k_ERK.*ERK.*Rosi.*Hs(ERK,ERK0_EEi,lambda_EEi,n_EEi);
df(15,:)=g_CEBP.*Rosi.*Hs(MEK,MEK_MC,lambda_MC,n_MC).*Hs(PPAR,PPAR_PC,lambda_PC,n_PC)-k_CEBP.*Hs(ERK,ERK0_EC,lambda_EC,n_EC).*CEBP.*TGF;
df(16,:)=g_PPAR.*Rosi.*Hs(MEK,MEK_MP,lambda_MP,n_MP).*Hs(CEBP,CEBP0_CP,lambda_CP,n_CP)-k_PPAR.*Hs(ERK,ERK0_EP,lambda_EP,n_EP).*PPAR;

function H=Hs(X,X0,lamda,n)
H=(1+lamda.*(X./X0).^n)./(1+(X./X0).^n);

