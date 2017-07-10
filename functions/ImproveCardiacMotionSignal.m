function para=ImproveCardiacMotionSignal(Cardiac_Signal,para);

%%% find the End diastole and End systole
[valley_values,valley_index]= findpeaks(-double(Cardiac_Signal));
valley_values=-valley_values;
%time = para.TR:para.TR:nt*2*para.TR;
plot(Cardiac_Signal,'r','LineWidth',1)
hold on

plot(valley_index,valley_values,'g*')
hold off

[peak_values,peak_index]= findpeaks(double(Cardiac_Signal));

index=0;
span=6;
while length(valley_index)>ceil(para.HeartFS*para.TR*para.nt)*1.2
    index=1;
    span=span+2;
    Cardiac_Signal_new = smooth(Cardiac_Signal,span,'lowess');
    [valley_values,valley_index]= findpeaks(-double(Cardiac_Signal_new));
    %[peak_values,peak_index]= findpeaks(double(Cardiac_Signal));
    valley_values=-valley_values;
    plot(Cardiac_Signal_new,'r','LineWidth',1)
    hold on
    plot(valley_index,valley_values,'g*')
    hold off
end
if index==1
    Cardiac_Signal=Cardiac_Signal_new; clear Cardiac_Signal_new;
    [valley_values,valley_index]= findpeaks(-double(Cardiac_Signal));
    [peak_values,peak_index]= findpeaks(double(Cardiac_Signal));
    valley_values=-valley_values;
end
clear index
plot(Cardiac_Signal,'r','LineWidth',1)
hold on
plot(valley_index,valley_values,'g*')
hold off

%%%exclude local false minima 
for ii=1:length(valley_index)
    t1=valley_values(ii)<Cardiac_Signal(valley_index(ii)+1 : min(valley_index(ii)+ceil(1/(para.HeartFS*para.TR)*0.2),length(Cardiac_Signal)));
    t2=valley_values(ii)<Cardiac_Signal(max(valley_index(ii)-ceil(1/(para.HeartFS*para.TR)*0.2),1) : valley_index(ii)-1);
    if isempty(find([t1;t2]==0));
        valley_index_new(ii)=valley_index(ii);
    else
        valley_index_new(ii)=0;
    end
end

%%%exclude local false maxima 
for ii=1:length(peak_index)
    t1=peak_values(ii)>Cardiac_Signal(peak_index(ii)+1 : min(peak_index(ii)+ceil(1/(para.HeartFS*para.TR)*0.2),length(Cardiac_Signal)));
    t2=peak_values(ii)>Cardiac_Signal(max(peak_index(ii)-ceil(1/(para.HeartFS*para.TR)*0.2),1) : peak_index(ii)-1);
    if isempty(find([t1;t2]==0));
        peak_index_new(ii) = peak_index(ii);
    else
        peak_index_new(ii) = 0;
    end
end

ES_index=valley_index(find(valley_index_new~=0));
ES_values=valley_values(find(valley_index_new~=0));
ED_index=peak_index(find(peak_index_new~=0));
ED_values=peak_values(find(peak_index_new~=0));
clear valley_index_new peak_index_new
plot(Cardiac_Signal,'r','LineWidth',1)
hold on
plot(ES_index,ES_values,'g*')
hold off

para.CardiacPhase=ceil(mean(diff(ES_index)));
para.Sys=ceil(para.CardiacPhase*0.4); %%%% Change to 0.5 for very long cycles
para.Dia=para.CardiacPhase-para.Sys;

disp(sprintf('para.CardiacPhase: %f', para.CardiacPhase));
disp(sprintf('para.Sys: %f', para.Sys));

if ES_index(1)<para.Sys;
    ES_index=ES_index(2:end);
    ES_values=ES_values(2:end);
end
if ES_index(end)+para.Dia+1>para.nt
    ES_index=ES_index(1:end-1);
    ES_values=ES_values(1:end-1);
end
while ED_index(2)< ES_index(1)
    ED_index=ED_index(2:end);
    ED_values=ED_values(2:end);
end

id = 1;
%find average sys:
for ii=1:length(valley_index)
  while (peak_index(id) < valley_index(ii)) && (id < length(peak_index))
      id = id + 1;
  end
  if (id > 1)
    sys_length(ii) = valley_index(ii) - peak_index(id-1);
  else
    sys_length(ii) = 0;
  end
end

sys_length1 = sys_length(find(sys_length~=0))
disp(sprintf('other sys: %f', mean(sys_length1)));
sys_length1


plot(Cardiac_Signal,'r','LineWidth',1)
hold on
plot(ES_index,ES_values,'g*')
plot(ED_index,ED_values,'b+')
hold off
para.ES_index=ES_index;
para.ED_index=ED_index;
para.ntres=length(para.ES_index);