function para=ImproveCardiacMotionSignal(Cardiac_Signal,para);

%%% find the End diastole and End systole
[valley_values,valley_index]= findpeaks(-double(Cardiac_Signal));
valley_values=-valley_values;
%time = para.TR:para.TR:nt*2*para.TR;
plot(Cardiac_Signal,'r','LineWidth',1)
hold on

plot(valley_index,valley_values,'g*')
hold off

index=0;
span=para.span;
while length(valley_index)>ceil(para.HeartFS*para.TR*para.nt)+4
    index=1;
    span=span+2;
    Cardiac_Signal_new = smooth(Cardiac_Signal,span,'lowess');
    [valley_values,valley_index]= findpeaks(-double(Cardiac_Signal_new));
    valley_values=-valley_values;
    plot(Cardiac_Signal_new,'r','LineWidth',1)
    hold on
    plot(valley_index,valley_values,'g*')
    hold off
end
if index==1
    Cardiac_Signal=Cardiac_Signal_new; clear Cardiac_Signal_new
end
clear index
[valley_values,valley_index]= findpeaks(-double(Cardiac_Signal));
valley_values=-valley_values;
plot(Cardiac_Signal,'r','LineWidth',1)
hold on
plot(valley_index,valley_values,'g*')
hold off

%%%exclude local false minimal 
for ii=1:length(valley_index)
    t1=valley_values(ii)<Cardiac_Signal(valley_index(ii)+1:min(valley_index(ii)+ceil(1/(para.HeartFS*para.TR)/2),para.nt));
    t2=valley_values(ii)<Cardiac_Signal(max(valley_index(ii)-ceil(1/(para.HeartFS*para.TR)/2),1):valley_index(ii)-1);
    if isempty(find([t1;t2]==0));
        valley_index_new(ii)=valley_index(ii);
    else
        valley_index_new(ii)=0;
    end
end

ES_index=valley_index(find(valley_index_new~=0));
ES_values=valley_values(find(valley_index_new~=0));
clear valley_index_new peak_index_new
plot(Cardiac_Signal,'r','LineWidth',1)
hold on
plot(ES_index,ES_values,'g*')
hold off

para.CardiacPhase=ceil(mean(diff(ES_index)));
para.Sys=ceil(para.CardiacPhase*0.4);
para.Dia=para.CardiacPhase-para.Sys;

if ES_index(1)<para.Sys;
    ES_index=ES_index(2:end);
    ES_values=ES_values(2:end);
end
if ES_index(end)+para.Dia+1>para.nt
    ES_index=ES_index(1:end-1);
    ES_values=ES_values(1:end-1);
end
plot(Cardiac_Signal,'r','LineWidth',1)
hold on
plot(ES_index,ES_values,'g*')
hold off
para.ES_index=ES_index;
para.ntres=length(para.ES_index);