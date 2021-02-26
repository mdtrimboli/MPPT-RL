Limit_V=25;
Limit_P=70;
T = 10:1:45;
G = 1;
x_voc = -5:0.1:Limit_V;
count=1;

Ph=zeros(size(T,2)*size(x_voc,2),1);
Vh=zeros(size(T,2)*size(x_voc,2),1);
Gh=zeros(size(T,2)*size(x_voc,2),1);


for k=1:size(T,2)
    fprintf('k=%i\n',k);
    for i=1:size(x_voc,2)
        [P,V]=altpvmodel(G,T(k),x_voc(i));
        Ph(count,1)=P;
        Vh(count,1)=V;
        Gh(count,1)=G;
        count=count+1;
    end
end

Vh=round(Vh*10);
state_list2=[Gh Vh Ph];
state_list = unique(state_list2,'rows');

%Generate a state list
% states=zeros(length(x1)*length(x2)*length(x3),2); % 2 Column matrix of all possible combinations of the discretized state.
% index=1;
% for j=1:length(x1)
%     for k = 1:length(x2)
%         for l = 1:length(x3)
%         states(index,1)=x1(j);
%         states(index,2)=x2(k);
%         states(index,3)=x3(l);
%         index=index+1;
%         end
%     end
% end
