clear;close all;


% RITHIKA CHANGE ME %
ecgSamplingRate = 977;
beginTime_msec = -45; % time before R peak to cut 
endTime_msec = 575; % time after R peak to cut 

path = '/home/dan/Desktop/rithika waterfall plots/Raw_ECG_V_timeseries';
filename = 'X001B.bin';

keepIdxs = 

% [filename, path] = uigetfile('*.txt;*.bin');
channel = 2;
if contains(filename,'.bin')
    disp(['Loading ' filename ' channel ' num2str(channel)]);
    fileID = fopen(fullfile(path,filename));
    data = fread( fileID, [ 23  , 517154*23 ], 'double')' ;
elseif contains(filename, '.txt')
    disp(['Loading ' filename ' channel ' num2str(channel)]);
    data = readmatrix(fullfile(path,filename), 'Delimiter', ' ');
else
    disp(['File ' filename ' not suported. Do .txt or .bin']);
end
ecg = data(:,channel);
output_name = filename(1:end-4); % will append a _processed and write out mat/csv of arrays



%ecg = smooth(ecg,11,'sgolay');
ecg(isnan(ecg))=nanmean(ecg);

rs = DeepR_and_PVC(ecg,ecgSamplingRate);
% find the peak around where the r wave is. Use raw ecg now

rpeaks = RPeakfromRWave(ecg,rs);

% plot(ecg);
% hold on;
% plot(rpeaks,ecg(rpeaks), '*');

% convert msec to samples
beginTime_msec = beginTime_msec * (ecgSamplingRate/1000);
endTime_msec = endTime_msec * (ecgSamplingRate/1000);

% create begining points to start beat
begining_cutoff = round(rpeaks + beginTime_msec);
begining_cutoff_mask = begining_cutoff < 0; % we don't want anything remove any <= 0. create a array where true if < 0

% create end points to end beat
end_cutoff = round(rpeaks + endTime_msec);
end_cutoff_mask = end_cutoff > length(ecg); % we don't want anything remove any > length. create a array where true if > length

% remove any intervals > length or < 0
begining_cutoff(begining_cutoff_mask | end_cutoff_mask) = [];
end_cutoff(begining_cutoff_mask | end_cutoff_mask) = [];


% calculate longest beat
longest = max(end_cutoff-begining_cutoff);

% create an array that has rows == number of beats and columns == longest
rawarray = zeros(length(end_cutoff),round(longest));

% fill in array. Subtract min point before r peak to kinda "align" amplitudes by q.
for x = 1:length(end_cutoff)
    beat = ecg(begining_cutoff(x):end_cutoff(x));
    [~,idx] = max(beat);
    [~,minidx] = min(beat(1:idx));
    beat = beat - beat(minidx);
    rawarray(x,1:length(beat)) = beat;
end


%% plot all pretty ugly because of outliers
figure();
% plot(rawarray');
hold on;
for x = 1:size(rawarray,1)
    red = plot3(1:size(rawarray,2),rawarray(x,:),ones(size(rawarray,2),1)*x, 'r');
    while true
        x = input('Keep?');
        if strcmp(x,'y')
            set(red, 'Color', 'b')
            break
        elseif strcmp(x,'n')
            delete(red)
            break
        end
    end
end
title('All Data w/ outliers');


% plot outlierless in both amplitude/length. 
lenOutliers = ~isoutlier(end_cutoff-begining_cutoff,'mean'); % this makes an array of 0's/1's. 1's represent outliers.  The ~ in front inverts it so we get the 'good' ones

% this applys a function to each element of the array passed in.  In this
% case I pass in an array that is from 1 to the length of the array (ie
% 1:length(rawarray)).  The 'function' is @(ROWIDX) max(array(ROWIDX,:)) think
% of the @ part as your input variables and the next portion as what
% happens. so each element that is passed in from the initial array is
% 1,2,3,... then it takes the max of the row for the matrix that contains
% the ecg data.
maxAmp = arrayfun(@(ROWIDX) max(rawarray(ROWIDX,:)), 1:size(rawarray,1)); 
ampOutliers = ~isoutlier(maxAmp,'mean'); % this makes an array of 0's/1's. 1's represent outliers.  The ~ in front inverts it so we get the 'good' ones

%% plot all pretty ugly because we found a lot of beats
% figure();
% plot(rawarray(lenOutliers & ampOutliers,:)'); % use the 0/1 arrays then the and operater to get only the places where both arrays have 1's
% xlim([1,endTime_msec]);
% title('All Data w/out amplitude/length outliers');


% smoothed = rawarray;
% for x = 1:size(rawarray,1)
%     smoothed(x,:) = smooth(smoothed(x,:),11,'sgolay');
% end

%% plot "3d" waterfall without outliers
figure();
stackedplot(rawarray',3,1,4);
ylim([1,endTime_msec]);
zlim([-0.5, 2.0]);
xlim([1,50]);
alpha 0.0
ylabel('Time (msec)')
xlabel('Beat #')
zlabel('Voltage (mV)')
view(-90,0)
set(gca, 'fontsize', 18);
set(gca, 'fontweight', 'bold');
mymap = [ 0.000, 0.000, 0.000; 0.000, 0.000, 0.294; 0.000, 0.000, 0.588; 0.000, 0.000, 0.794; 0.000, 0.000, 1.000; 0.000, 0.098, 0.892; 0.000, 0.196, 0.784; 0.000, 0.294, 0.686; 0.000, 0.392, 0.588; 0.000, 0.588, 0.588; 0.000, 0.784, 0.588; 0.196, 0.784, 0.294; 0.392, 0.784, 0.000; 0.588, 0.784, 0.000; 0.784, 0.784, 0.000; 0.784, 0.539, 0.000; 0.784, 0.294, 0.000; 0.892, 0.147, 0.098; 1.000, 0.000, 0.196; 0.941, 0.000, 0.392; 0.882, 0.000, 0.588; 0.784, 0.000, 0.686; 0.686, 0.000, 0.784; 0.637, 0.000, 0.838; 0.588, 0.000, 0.892; 0.539, 0.000, 0.946; 0.490, 0.000, 1.000 ];
colormap(mymap)
%caxis([-0.5, 1.8])



