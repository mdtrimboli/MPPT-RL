%% MPPT Comparison
%SARSA vs QL vs PyO vs Incremental Conductance

Limit_V=25;
x_voc = 0:0.01:Limit_V;
actions = [-3 -0.5 -0.1 0.1 0.5 3];
action=0.25;
iter_max=50000;
step_max=35;
count=0;
T = 10:1:45;
G = 0.1:0.01:1;

    V_hist=zeros(4,iter_max*step_max); % F1=P&O, F2=IC, F3=QL, F4=S
    P_hist=zeros(4,iter_max*step_max); % F1=P&O, F2=IC, F3=QL, F4=S
    MPP_hist=zeros(1,iter_max*step_max);
    
    Error_ab=zeros(4,iter_max*step_max);
    Errorab_hist_5=zeros(4,iter_max*step_max);
    E_hist_per=zeros(4,iter_max*step_max);
    E_hist_ij=zeros(4,iter_max*step_max);
    Error_hist_5=zeros(4,iter_max);
    action_hist_QL=zeros(1,iter_max*step_max);
    action_hist_S=zeros(1,iter_max*step_max);
    rastreo=zeros(1,step_max);
%% Load States Space
load('state_list_e3.mat');
%% Initials Conditions for QL
load('QT_QL_e3.mat')
%% Initials Conditions for SARSA
load('QT_S_e3.mat')
%% External Loop

sel_G=randi(length(G),iter_max,1);
sel_T=randi(length(T),iter_max,1);
e=rand(1,iter_max);   


for j=1:iter_max
    fprintf('epoch=%i\n',j);
    Go=G(sel_G(j));
    To=T(sel_T(j));
    [MPP, Vmpp]=MPPT(To,Go);
    [Ph,Vh,Ih]=altpvmodel(Go,To,x_voc);
    Possible_Voc = fliplr(find(Ih>=0));
    Voc=x_voc(Possible_Voc(1));
    V=Voc*rand(1);
    Vr=round(V,1);
    [Po,Vo,Io]=altpvmodel(Go,To,Vr); %(G,T,V)

    Vo_Pyo=Vo;
    Po_Pyo=Po;
    delta_P_Pyo=0;
    delta_V_Pyo=0;
    
    Vo_IC=Vo;
    Io_IC=Io;
    Po_IC=Po;
    delta_V_IC=0;
    delta_I_IC=0.01;
    Pend=0;
    
    Vo_QL=Vo;
    Po_QL=Po;
    
    Vo_S=Vo;
    Po_S=Po;
    
     
 for i=1:step_max
    %% Internal Loop
 
    count=count+1;  
    [Vn_Pyo,Pn_Pyo,delta_P_Pyo,delta_V_Pyo] = PyO(Po_Pyo,Vo_Pyo,delta_P_Pyo,delta_V_Pyo,action,G(sel_G(j)),T(sel_T(j)));
    Vo_Pyo=Vn_Pyo;
    Po_Pyo=Pn_Pyo;
    
    [Vn_IC,Pn_IC,In_IC,Pend,delta_V_IC,delta_I_IC] = IC_control(Vo_IC,Io_IC,Po_IC,Pend,delta_V_IC,delta_I_IC,action,G(sel_G(j)),T(sel_T(j)));
    Vo_IC=Vn_IC;
    Po_IC=Pn_IC;
    Io_IC=In_IC;

     [Vn_QL,Pn_QL,act_QL]=test_QL_e3(Vo_QL,Po_QL,state_list,actions,Q_table,Go,To);
     action_hist_QL(count)=act_QL;
     Vo_QL = Vn_QL;
     Po_QL = Pn_QL;
    
     [Vn_S,Pn_S,act_S]=test_Sarsa_e3(Vo_S,Po_S,state_list,actions,Sarsa_table,G(sel_G(j)),T(sel_T(j)));
     action_hist_S(count)=act_S;
     Vo_S = Vn_S;
     Po_S = Pn_S;

    %% Captura de datos por iteracion
    V_hist(:,count)=[Vo_Pyo; Vo_IC; Vo_QL; Vo_S];
    P_hist(:,count)=[Po_Pyo; Po_IC; Po_QL; Po_S];
    
    MPP_hist(count)=MPP;
 
    Error_ab(:,count)=[abs(MPP-Pn_QL); abs(MPP-Pn_S); abs(MPP-Pn_Pyo); abs(MPP-Pn_IC)];
    E_hist_ij(:,count)=[MPP-Pn_QL; MPP-Pn_S; MPP-Pn_Pyo; MPP-Pn_IC];
    E_hist_per(:,count)=[sqrt((E_hist_ij(1,count).^2))*100./MPP; sqrt((E_hist_ij(2,count).^2))*100./MPP; ...
        sqrt((E_hist_ij(3,count).^2))*100./MPP; sqrt((E_hist_ij(4,count).^2))*100./MPP];
    
 end
 
%     scatter(Vo_Pyo,Po_Pyo,'filled','green')
%     hold on
%     scatter(Vo_IC,Po_IC,'filled','blue')
%     scatter(Vo_QL,Po_QL,'filled','yellow')
%     scatter(Vo_S,Po_S,'filled','red')
%     plot(x_voc,Ph,'black')
%     scatter(Vmpp,MPP,'black')
%     hold off
%     axis([0 25 0 70])
%     xlabel('Voltage')
%     ylabel('Power')
 
end

MAE_QL_t = sum(Error_ab(1,:))./(step_max*iter_max);
MAE_S_t = sum(Error_ab(2,:))./(step_max*iter_max);
MAE_PO_t = sum(Error_ab(3,:))./(step_max*iter_max);
MAE_IC_t = sum(Error_ab(4,:))./(step_max*iter_max);

MSE_QL_t = sum(E_hist_ij(1,:).^2)./(step_max*iter_max);
MSE_S_t = sum(E_hist_ij(2,:).^2)./(step_max*iter_max);
MSE_PO_t = sum(E_hist_ij(3,:).^2)./(step_max*iter_max);
MSE_IC_t = sum(E_hist_ij(4,:).^2)./(step_max*iter_max);

fprintf('MSE_QL=%f\n',MSE_QL_t);
fprintf('MSE_S=%f\n',MSE_S_t);
fprintf('MSE_PO=%f\n',MSE_PO_t);
fprintf('MSE_IC=%f\n\n',MSE_IC_t);

fprintf('MAE_QL=%f\n',MAE_QL_t);
fprintf('MAE_S=%f\n',MAE_S_t);
fprintf('MAE_PO=%f\n',MAE_PO_t);
fprintf('MAE_IC=%f\n\n',MAE_IC_t);


c5=0;
for k=1:(step_max*iter_max)
    if mod(k,i)==0
        c5=c5+1;
        Error_hist_5(1,c5)=E_hist_ij(1,k);
        Error_hist_5(2,c5)=E_hist_ij(2,k);
        Error_hist_5(3,c5)=E_hist_ij(3,k);
        Error_hist_5(4,c5)=E_hist_ij(4,k);
        
        Errorab_hist_5(1,c5)=Error_ab(1,k);
        Errorab_hist_5(2,c5)=Error_ab(2,k);
        Errorab_hist_5(3,c5)=Error_ab(3,k);
        Errorab_hist_5(4,c5)=Error_ab(4,k);
    end
end
MSE_QL_pe = sum(Error_hist_5(1,:).^2)/j;
MSE_S_pe= sum(Error_hist_5(2,:).^2)/j;
MSE_PO_pe = sum(Error_hist_5(3,:).^2)/j;
MSE_IC_pe = sum(Error_hist_5(4,:).^2)/j;

MAE_QL_pe = sum(Errorab_hist_5(1,:))/j;
MAE_S_pe= sum(Errorab_hist_5(2,:))/j;
MAE_PO_pe = sum(Errorab_hist_5(3,:))/j;
MAE_IC_pe= sum(Errorab_hist_5(4,:))/j;

fprintf('MSE-QL (Fin episodio)=%f\n',MSE_QL_pe);
fprintf('MSE-SARSA (Fin episodio)=%f\n',MSE_S_pe);
fprintf('MSE-P&O (Fin episodio)=%f\n',MSE_PO_pe);
fprintf('MSE-IC (Fin episodio)=%f\n\n',MSE_IC_pe);

fprintf('MAE-QL (Fin episodio)=%f\n',MAE_QL_pe);
fprintf('MAE-SARSA (Fin episodio)=%f\n',MAE_S_pe);
fprintf('MAE-P&O (Fin episodio)=%f\n',MAE_PO_pe);
fprintf('MAE-IC (Fin episodio)=%f\n\n',MAE_IC_pe);