%% MPPT using Q-Learning
% Author: Maximiliano Trimboli (FICA-UNSL)
%% Define variables
Limit_sup_V=25;
Limit_inf_V=-5;
Value=0;
Cont_error=0;
Ph=zeros(1,Limit_sup_V);
Vh=zeros(1,Limit_sup_V);
Ih=zeros(1,Limit_sup_V);
x_voc = 0:0.1:Limit_sup_V;
%% State space
load('state_list_e3.mat');
states=state_list;
state_list=0;

actions = [-3 -0.5 -0.1 0.1 0.5 3];
T = 10:1:45;
G = 0.1:0.01:1;
s_states = size(states);
size_states = s_states(1);
s_actions = size(actions);
size_actions = s_actions(2);

%% Define Q table and reward
Q_table = zeros(size_states, size_actions);

%% Initialize parameters
iter_max = 500000;
steps=35;
gamma = 0.5; %Discount factor [Default 0.9]
alpha = 0.5; %Initial learning rate
epsilon = 1; % Epsilon greedy [Default 0.9]
epsilon_decay = 0.999999; %Decay factor per iteration

Voltage_hist=zeros(1,iter_max);
Power_hist=zeros(1,iter_max);
reward_s=zeros(1,10);
E_hist=zeros(1,iter_max);
EA_hist=zeros(1,iter_max);
ERP_hist=zeros(1,iter_max);
Average=zeros(1,iter_max);
V_hist=zeros(1,iter_max);
MPP_hist=zeros(1,iter_max);
%% QL Algorithm

sel_T=randi(length(T),iter_max,1);
sel_G=randi(length(G),iter_max,1);
e=rand(1,iter_max);  

for i=1:iter_max
     fprintf('Epoch %i\n',i); 
    %Lectura de datos de variables P,V,dV [Armar instrucciones para eso]
    %Initial conditions of simulation
    Go=G(sel_G(i));
    To=T(sel_T(i));
    MPP=MPPT(T(sel_T(i)),G(sel_G(i)));
    [Ph,Vh,Ih]=altpvmodel(G(sel_G(i)),T(sel_T(i)),x_voc);
    Possible_Voc = fliplr(find(Ih>=0));
    Voc=x_voc(Possible_Voc(1));
    V=Voc*rand(1);
    Vr=round(V,1);
    [Po,Vo,Io]=altpvmodel(G(sel_G(i)),T(sel_T(i)),Vr); %(G,T,V)
    
    % Search initial state
    Vr_n=int16(Vr*10);
    index_ist = find(states(:,1)==Go); %List of indexes of the possible next states
    selected_states_t = [index_ist, states(index_ist,2), states(index_ist,3)];
    index_ns = find(selected_states_t(:,2)==Vr_n); 
    selected_states = [index_ns, selected_states_t(index_ns,3)];
    [~,s_index2]=min(abs(selected_states(:,2)-Po));
    s_index = selected_states(s_index2,1);
    
    j=0;
    Value=0;
  
    for step=1:steps

        %Epsilon greedy policy       
        if e(i)<epsilon
            %Exploration
            a_index=randi([1 length(actions)],1,1);
            act_sel=actions(a_index);
%             fprintf('Exploration (epsilon:%.3f)\n',epsilon);
        else
            %Exploitation
            [max_q,a_index] = max(Q_table(s_index,:));
            max_index=find(Q_table(s_index,:)==max_q);
            a_index=max_index(randi(length(max_index)));
            act_sel=actions(a_index);    
%             fprintf('Exploitation (epsilon:%.3f)\n',epsilon);
        end
    
        Vf = Vo + act_sel;
        [Pf,Vf,If] = altpvmodel(G(sel_G(i)),T(sel_T(i)),Vf);
        dV=Vf-Vo;
        dP=Pf-Po;
  
    %Reward
     %Rewards utilizados: 1)Reward fijo, 2)Funci�n de Pf, 3)Funci�n de dP
        if Pf>0
            reward=2*dP;
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
    
    
    %Define state index
        Vf_m=int16(Vf*10);
        index_ns = find(selected_states_t(:,2)==Vf_m); %List of indexes of the possible next states
        selected_states = [index_ns, selected_states_t(index_ns,3)];
        [~,index_ns2]=min(abs(selected_states(:,2)-Pf));
        next_state = selected_states(index_ns2,1);
    
    % Update q-values as per Bellman's equation
        maxQ=max(Q_table(next_state,a_index));
        Q_table(s_index,a_index)= Q_table(s_index,a_index)+alpha*(reward+gamma*maxQ-Q_table(s_index,a_index));
        epsilon=epsilon*epsilon_decay;
        s_index=next_state;
        
        Po=Pf;
        Vo=Vf;
    end
        Voltage_hist(i)=Vf;
        Power_hist(i)=Pf;
     if i<11
         reward_s(i)= Value;
    else
        reward_s(1)= reward_s(2);
        reward_s(2)= reward_s(3);
        reward_s(3)= reward_s(4);
        reward_s(4)= reward_s(5);
        reward_s(5)= reward_s(6);
        reward_s(6)= reward_s(7);
        reward_s(7)= reward_s(8);
        reward_s(8)= reward_s(9);
        reward_s(9)= reward_s(10);
        reward_s(10)= Value;
    end
    
    Error=(MPP-Pf);
    Error_abs=abs(Error);
    E_hist(i)=Error;
    EA_hist(i)=Error_abs;
    ERP_hist(i)=Error_abs.*100./MPP;
    
    size_rew = size(reward_s);
    Average_s = sum(reward_s)/size_rew(2);
    Average(i)= Average_s;
    V_hist(i)=Value;  
end

%% Metrics
Mean_sq_error = sum(E_hist.^2)/iter_max;
fprintf('MSE=%.2f\n\n\n',Mean_sq_error);
%% Plots
c = linspace(1,10,length(Voltage_hist));

figure
plot(1:iter_max,Average)
xlabel('Episode Number')
ylabel('Reward')
% legend(alpha_st)

figure
scatter(Voltage_hist,Power_hist, [],c)
xlabel('Voltage')
ylabel('Power')


