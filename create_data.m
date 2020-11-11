clear;close all;

% setup variables

filename = '/home/dan/git/work/waterfallPlot/HPLN-R14-2 3-11-19_JC.txt';
filename = '/home/dan/Downloads/Raw_ECG_V_timeseries/X000.txt';
output_name = 'rithika'; % will append a _processed and write out mat/csv of arrays
ecgSamplingRate = 977;


beginTime_msec = -45; % time before R peak to cut
endTime_msec = 600; % time after R peak to cut


% read in tabbed file
data = readmatrix(filename, 'Delimiter', ' ');

% ecg == 2nd column
for i=1:1
    ecg = data(:,2);
    ecg(ecg>1) = 1;
    ecg(ecg<-.4) = -1;
    ecg_outlierless = filloutliers(ecg,'nearest','mean');
    rs = rpeakdetect(ecg_outlierless,ecgSamplingRate);

    % find the peak around where the r wave is. Use raw ecg now
    rpeaks = RPeakfromRWave(ecg,rs);
    figure();

    plot(ecg);
    hold on;
    plot(rpeaks,ecg(rpeaks),'or');
end

% make nans == mean of wave
disp(['Nan count:' num2str(sum(isnan(ecg)))]);
disp('Making nans == mean of wave');
ecg(isnan(ecg))=nanmean(ecg);

% simple outlier removal
ecg_outlierless = filloutliers(ecg,'nearest','mean');

% move ecg mean to ~ 0
ecg_outlierless = ecg_outlierless - mean(ecg_outlierless);

% simple r wave detector, not the best at discriminating noise, use
% modified ecg to detect
rs = rpeakdetect(ecg_outlierless,ecgSamplingRate);

% find the peak around where the r wave is. Use raw ecg now
rpeaks = RPeakfromRWave(ecg,rs);

% plot
% figure();
% plot(ecg); 
% hold on; 
% plot(rpeaks,ecg(rpeaks),'o');
% title('Raw ecg with detected r points');
% 
% figure();
% plot(ecg_outlierless); 
% hold on; 
% plot(rpeaks,ecg_outlierless(rpeaks),'o');
% title('Outlierless ecg with r points');


% convert msec to samples
beginTime_msec = beginTime_msec * (ecgSamplingRate/1000);
endTime_msec = endTime_msec * (ecgSamplingRate/1000);

% create begining points to start beat
begining_cutoff = rpeaks + beginTime_msec;
begining_cutoff_mask = begining_cutoff < 0; % we don't want anything remove any <= 0. create a array where true if < 0

% create end points to end beat
end_cutoff = rpeaks + endTime_msec;
end_cutoff_mask = end_cutoff > length(ecg); % we don't want anything remove any > length. create a array where true if > length

% remove any intervals > length or < 0
begining_cutoff(begining_cutoff_mask | end_cutoff_mask) = [];
end_cutoff(begining_cutoff_mask | end_cutoff_mask) = [];

% plot first few 'beats' ~ 2 sec and begin (green) to end (red). Used this
% to get -20 and 100 numbers in very beginning
% figure();
% end_val = 2000;
% sec2_ecg = ecg(1:end_val);
% plot(sec2_ecg);
% hold on;
% sec2_ecg_rpeaks = rpeaks(rpeaks < end_val);
% plot(sec2_ecg_rpeaks,sec2_ecg(sec2_ecg_rpeaks),'o');
% vline(begining_cutoff(1:20),'g');
% vline(end_cutoff(1:20),'r');
% xlim([1,end_val]);

% calculate longest beat
longest = ceil(max(end_cutoff-begining_cutoff));

% create an array that has rows == number of beats and columns == longest
% beat
rawarray = zeros(length(end_cutoff),longest);

% fill in array. Subtract min point before r peak to kinda "align" amplitudes by q.
for x = 1:length(begining_cutoff)
    beat = ecg(begining_cutoff(x):end_cutoff(x));
    [~,idx] = max(beat);
    [~,minidx] = min(beat(1:idx));
    beat = beat - beat(minidx);
    rawarray(x,1:length(beat)) = beat;
end


%% plot all pretty ugly because of outliers
% figure();
% plot(rawarray');
% xlim([1,endTime_msec]);
% title('All Data w/ outliers');

% plot outlierless in both amplitude/length. 
lenOutliers = ~isoutlier(end_cutoff-begining_cutoff,'mean'); % this makes an array of 0's/1's. 1's represent outliers.  The ~ in front inverts it so we get the 'good' ones

% this applys a function to each element of the array passed in.  In this
% case I pass in an array that is from 1 to the length of the array (ie
% 1:length(rawarray)).  The 'function' is @(ROWIDX) max(array(ROWIDX,:)) think
% of the @ part as your input variables and the next portion as what
% happens. so each element that is passed in from the initial array is
% 1,2,3,... then it takes the max of the row for the matrix that contains
% the ecg data.
% maxAmp = arrayfun(@(ROWIDX) max(rawarray(ROWIDX,:))-min(rawarray(ROWIDX,:)), 1:length(rawarray)); 
% ampOutliers = ~isoutlier(maxAmp,'mean'); % this makes an array of 0's/1's. 1's represent outliers.  The ~ in front inverts it so we get the 'good' ones


compositeMean = NaN(1, size(rawarray, 2));
compositeSTD = NaN(1, size(rawarray, 2));

for I = 1:size(rawarray, 2)
    compositeMean(I) = nanmean(rawarray(:, I));
    compositeSTD(I) = nanstd(rawarray(:, I));
end
waveScore = zeros(1, size(rawarray, 1));

for B = 1:size(rawarray, 1)

    for I = 1:size(rawarray, 2)

        if ((rawarray(B, I) > (compositeMean(I) + 2.5 * compositeSTD(I)) || rawarray(B, I) < (compositeMean(I) - 2.5 * compositeSTD(I))))
            waveScore(B) = waveScore(B) + 1;
        end

        if ((rawarray(B, I) > (compositeMean(I) + 3 * compositeSTD(I)) || rawarray(B, I) < (compositeMean(I) - 3 * compositeSTD(I))))
            waveScore(B) = waveScore(B) + 3;
        end

        if ((rawarray(B, I) > (compositeMean(I) + 4 * compositeSTD(I)) || rawarray(B, I) < (compositeMean(I) - 4 * compositeSTD(I))))
            waveScore(B) = waveScore(B) + 9;
        end

    end

end

waveScore = waveScore ./ size(rawarray, 2);
score = waveScore > .30;



%% plot all pretty ugly because we found a lot of beats
% figure();
% plot(rawarray(lenOutliers & ampOutliers,:)'); % use the 0/1 arrays then the and operater to get only the places where both arrays have 1's
% xlim([1,endTime_msec]);
% title('All Data w/out amplitude/length outliers');


%% plot "3d" waterfall without outliers
figure();

stackedplot(rawarray',3,10);
ylim([1,endTime_msec]);
alpha 0.0;
title('all');
% 
% figure();
% %mask = ~(lenOutliers & ampOutliers);
% mask = score;
% outlierless_array = rawarray;
% outlierless_array(mask,:) = 0;
% stackedplot(outlierless_array',3,10);
% alpha 0.0;
% title('zero for outlier');


% save([output_name '_processed.mat'], 'rawarray', 'outlierless_array');
% csvwrite([output_name '_allwaves.txt'],'rawarray');
% csvwrite([output_name '_outlierless.txt'],'outlierless_array');


