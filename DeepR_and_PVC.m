% try deep PVC envelope

function [HIll_climbedRpeaks]  = DeepR_and_PVC(Epoch_orig,sample_rate)
dataset = 'MESA';
dbg =  0;

% R peak identification for different data set
if strcmp('MESA',dataset) == 1
    mean_epoch = Epoch_orig - mean(Epoch_orig);
    Rpeaks = ecgLms(Epoch_orig - mean(Epoch_orig),sample_rate,5,12);
%     Rpeaks_other = ecgLms_allpeaks(Epoch_orig - mean(Epoch_orig),sample_rate,5,12);
elseif strcmp('SHHS1',dataset) == 1
    mean_epoch = Epoch_orig - mean(Epoch_orig);
    Rpeaks = ecgLms(mean_epoch,sample_rate,5,12,0);
%     Rpeaks_other = ecgLms_allpeaks(mean_epoch,sample_rate,5,12);
elseif strcmp('GP',dataset) == 1
    Rpeaks = ecgLmsDog1(onechanPLUS7 - mean(onechanPLUS7),sample_rate/2,5,12,'y');
end

HIll_climbedRpeaks = hill_climber(Epoch_orig,Rpeaks,sample_rate);

% % PVC detection
% VERYTEMPlistcycok = ecgSL(Rpeaks);
% tempgood_i = find(VERYTEMPlistcycok==1); 
% LLLL = diff(Rpeaks);   
% upto88 = min([ 20 numel(LLLL) ]);  
% if upto88<1
%    upto88=1; 
% end 
% startingAvgUnitsIndexes = median(LLLL(1:upto88));  
% prevStdMsec = 0 ;
% if numel(tempgood_i)>=15
%    tempgood_i = tempgood_i(1:15);
%    startingAvgUnitsIndexes = median(LLLL(tempgood_i)); 
%    prevStdMsec = std(LLLL(tempgood_i)*1000/125);
% end ; clear LLLL VERYTEMPlistcycok tempgood_i upto88 ;
% 
% LLLL = diff(Rpeaks);
% listcycok = ecgSLandEnvelope(Rpeaks,sample_rate,[10000 20000],startingAvgUnitsIndexes,prevStdMsec );
% 
% if dbg ==1
%     time = (1:length(Epoch_orig))/sample_rate;
%     figure;plot(time,Epoch_orig);hold on;plot(Rpeaks/sample_rate,Epoch_orig(Rpeaks),'o');
%     xlabel('time (s)');ylabel('ECG (mV)');
%     title('Original ECG with Barrys code R peaks'); 
% end

end