%% MPPT using SARSA
% Author: Maximiliano Trimboli (FICA-UNSL)
%% Define variables
Limit_sup_V=25;
Limit_inf_V=-5;
Value=0;
Ph=zeros(1,Limit_sup_V);
Vh=zeros(1,Limit_sup_V);
Ih=zeros(1,Limit_sup_V);
x_voc = 0:0.1:Limit_sup_V;

%% States space
load('state_list_e2.mat');
states=state_list;
state_list=0;

actions = [-3 -0.5 -0.1 0.1 0.5 3];
T = 10:1:45;
G = 1;
s_states = size(states);
size_states = s_states(1);
s_actions = size(actions);
size_actions = s_actions(2);

%% Define Q table and reward
Q_table = zeros(size_states, size_actions);

%% Initialize parameters
iter_max = 250000;
steps=35;
gamma = 0.2; %Discount factor [Default 0.9]
alpha = 0.5; %Initial learning rate
epsilon = 1; % Epsilon greedy [Default 0.9]
epsilon_decay = 0.999999; %Decay factor per iteration

Voltage_hist=zeros(1,iter_max);
Power_hist=zeros(1,iter_max);
E_hist=zeros(1,iter_max);
EA_hist=zeros(1,iter_max);
ERP_hist=zeros(1,iter_max);
Average=zeros(1,iter_max);
V_hist=zeros(1,iter_max);
MPP_hist=zeros(1,iter_max);
%% SARSA Algorithm

sel_G=1;
sel_T=randi(length(T),iter_max,1);
e=rand(1,iter_max);
Go=G(sel_G);

for i=1:iter_max
    fprintf('Epoch %i\n',i);
    %Lectura de datos de variables P,V,dV [Armar instrucciones para eso]
    %Initial conditions of simulation
    To=T(sel_T(i));
    MPP=MPPT(To,Go);
    [Ph,Vh,Ih]=altpvmodel(Go,To,x_voc);
    Possible_Voc = fliplr(find(Ih>=0));
    Voc=x_voc(Possible_Voc(1));
    V=Voc*rand(1);
    Vr=round(V,1);
    [Po,Vo,Io]=altpvmodel(Go,T(sel_T(i)),Vr); %(G,T,V)
   
   % Search initial state
    Vr_n=int16(Vr*10);
    index_ist = find(states(:,1)==Go); %List of indexes of the possible next states
    selected_states_t = [index_ist, states(index_ist,2), states(index_ist,3)];
    index_isv = find(selected_states_t(:,2)==Vr_n);
    selected_states_v = [index_isv, selected_states_t(index_isv,3)];
    [~,s_index2]=min(abs(selected_states_v(:,2)-Po));
    actual_state_index = selected_states_v(s_index2,1);
    
    j=0;
    Value=0;

   %Epsilon greedy policy
    if e(i)<epsilon
      %Exploration
      actual_act_index=randi([1 length(actions)],1,1);
      actual_action=actions(actual_act_index);
%       fprintf('Exploration\n');
    else
       %Exploitation
      [max_q,actual_act_index] = max(Q_table(actual_state_index,:));
      max_index=find(Q_table(actual_state_index,:)==max_q);
      actual_act_index=max_index(randi(length(max_index)));
      actual_action=actions(actual_act_index);    
%       fprintf('Exploitation\n');
    end
    
    
    
% Steps loop
    for step=1:1:steps
        
        Vf = Vo + actual_action;
        [Pf,Vf,If] = altpvmodel(Go,T(sel_T(i)),Vf);
        dP=Pf-Po;
        dV=Vf-Vo;
  
    %Reward
     %Rewards utilizados: 1)Reward fijo, 2)Función de Pf, 3)Función de dP
        if Pf>0
            reward=dP;
        elseif Pf<=0 && dP>0
            reward=-5;
        else
            reward=-50;
        end
        
    %Cumulative reward
        Value=Value+reward*gamma^j; 
        j=j+1;
     
       %Evaluate state limits
        if Vf>Limit_sup_V
           Vf=Limit_sup_V;
        end
        if Vf<Limit_inf_V
           Vf=Limit_inf_V;
        end
        
    %Define next state index
        Vf_m=int16(Vf*10);
        index_isv = find(selected_states_t(:,2)==Vf_m); %List of indexes of the possible next states
        selected_states_v = [index_isv, selected_states_t(index_isv,3)]; % [Index selected for Vf, Possibles Pf]
        [~,s_index2]=min(abs(selected_states_v(:,2)-Pf));
        next_state_index = selected_states_v(s_index2,1);

        %Epsilon greedy policy  
    if e(i)<epsilon
      %Exploration
      next_act_index=randi([1 length(actions)],1,1);
      next_action=actions(next_act_index);
    else
      %Exploitation
      [max_q,next_act_index] = max(Q_table(next_state_index,:));
      max_index=find(Q_table(next_state_index,:)==max_q);
      next_act_index=max_index(randi(length(max_index)));
      next_action=actions(next_act_index);    
    end
            
    % Update q-values as per Bellman's equation
        Q_table(actual_state_index,actual_act_index)= Q_table(actual_state_index,actual_act_index)+alpha*(reward+gamma*Q_table(next_state_index,next_act_index)-Q_table(actual_state_index,actual_act_index));    
        epsilon=epsilon*epsilon_decay; 
        actual_state_index=next_state_index;
        actual_action=next_action;
        actual_act_index=next_act_index;
        
        Po=Pf;
        Vo=Vf;
    end
        Voltage_hist(i)=Vf;
        Power_hist(i)=Pf;
 
    
    Error=(MPP-Pf);
    Error_abs=abs(Error);
    E_hist(i)=Error;
    EA_hist(i)=Error_abs;

    V_hist(i)=Value;  
end

%% Metrics
Mean_sq_error = sum(E_hist.^2)/iter_max;
fprintf('MSE=%.2f\n\n\n',Mean_sq_error);
%% Plots
c = linspace(1,10,length(Voltage_hist));

figure
scatter(Voltage_hist,Power_hist, [],c, '.')
xlabel('Voltage')
ylabel('Power')



