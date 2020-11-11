function [template, bad_segment_fig] = EnsembleTemplate(wave, segment_idxs, std_multiplier, type, dbg_plot, fig_name)

    segment_idxs = [segment_idxs(1:end - 1) segment_idxs(2:end)];
    segment_lengths = segment_idxs(:, 2) - segment_idxs(:, 1) + 1;
    %ensamble_classification = zeros(length(segment_idxs));

    %% bad lengths
    med_length = median(segment_lengths);
    length_min = round(med_length - std_multiplier * std(segment_lengths));
    length_max = round(med_length + std_multiplier * std(segment_lengths));
    badlengths = (segment_lengths < length_min) | (segment_lengths > length_max);  
    
    %% seperate segements from wave ignoring bad lengths
    good_segment_count = numel(badlengths(badlengths == 0));
    good_seg_lengths = segment_lengths(badlengths == 0);
    good_seg_ids = segment_idxs(badlengths == 0, :);
    good_segments = NaN(good_segment_count, max(good_seg_lengths) + 1);
    good_segment_peaks = NaN(good_segment_count, 2);
    
    for I = 1:size(good_segments, 1)
        good_segments(I, 1:good_seg_lengths(I)) = wave(good_seg_ids(I, 1):good_seg_ids(I, 2));
        [good_segment_peaks(I, 1), good_segment_peaks(I, 2)] = max(good_segments(I, :));
    end
    
    if dbg_plot == 1
        bad_type_count = numel(badlengths(badlengths == 1));
        bad_total_count = numel(badlengths(badlengths == 1));
        bad_segments_lengths = segment_lengths(badlengths == 1);
        bad_seg_ids = segment_idxs(badlengths == 1, :);
        badSegments = NaN(bad_total_count, max(bad_segments_lengths) + 1);

        for I = 1:bad_total_count
            badSegments(I, 1:bad_segments_lengths(I)) = wave(bad_seg_ids(I, 1):bad_seg_ids(I, 2));
        end
        bad_segment_fig = figure('Name', fig_name,'NumberTitle','off','visible','off');
        UpdateBadSegFig(bad_segment_fig, badSegments, rgb('Black'));
    end

    %% remove bad peaks by amplitude
    avepeak = mean(good_segment_peaks(:, 1));
    devpeak = std(good_segment_peaks(:, 1));
    badpeaks = ((good_segment_peaks(:, 1) > (avepeak + std_multiplier * 2 * devpeak)) | (good_segment_peaks(:, 1) < (avepeak - std_multiplier * 2 * devpeak)));

    if dbg_plot == 1
        bad_type_count = [bad_type_count;numel(badpeaks(badpeaks == 1))];
        bad_total_count = bad_total_count + numel(badpeaks(badpeaks == 1));
        badSegments = good_segments(badpeaks == 1, :);
        UpdateBadSegFig(bad_segment_fig, badSegments, rgb('Lime'));
    end
    
    good_seg_lengths = good_seg_lengths(badpeaks == 0);
    good_segments = good_segments(badpeaks == 0, :);
    good_segment_peaks = good_segment_peaks(badpeaks == 0, :);
        
    %% remove bad peaks by position
    avepeak = mean(good_segment_peaks(:, 2));
    devpeak = std(good_segment_peaks(:, 2));
    badpeaks = ((good_segment_peaks(:, 2) > (avepeak + std_multiplier * 3 * devpeak)) | (good_segment_peaks(:, 2) < (avepeak - std_multiplier * 3 * devpeak)));
    
    if dbg_plot == 1
        bad_type_count = [bad_type_count;numel(badpeaks(badpeaks == 1))];
        bad_total_count = bad_total_count + numel(badpeaks(badpeaks == 1));
        badSegments = good_segments(badpeaks == 1, :);
        UpdateBadSegFig(bad_segment_fig, badSegments, rgb('Khaki'));
    end
    
    good_segment_count = numel(badpeaks(badpeaks == 0));
    good_seg_lengths = good_seg_lengths(badpeaks == 0);
    good_segments = good_segments(badpeaks == 0, :);
    good_segment_peaks = good_segment_peaks(badpeaks == 0, :);


%     %% max slope based pruning
%     slope = zeros(1, size(good_segments, 1));
% 
%     for B = 1:size(good_segments, 1)
%         slope(B) = max(diff(good_segments(B, :)));
%     end
% 
%     badslope = (slope > (mean(slope) + 1 * std(slope)))';
% 
%     if dbg == 1
%         bad_total_count = bad_total_count + size(badslope(badslope == 1), 1);
%         badSegments = good_segments(badslope == 1, :);
%         UpdateBadSegFig(bad_segment_fig, badSegments, rgb('DarkSlateGray'));
%     end
%     
%     good_segment_count = size(badslope(badslope == 0), 1);
%     good_seg_lengths = good_seg_lengths(badslope == 0);
%     good_segments = good_segments(badslope == 0, :);
%     
%     figure();
%     plot(good_segments','b');
%     hold on;
%     plot(badSegments','r');

%     for i = 1:length(good_segments)
%         
%     end
    

    
    
    %% align
    if strcmpi(type,'ppg')
        good_seg_diff_peak = NaN(good_segment_count, 2);
        %% by foot
        %[good_seg_diff_peak(:,1),good_seg_diff_peak(:,2)] = find_foot_pulseox(good_segments,0);
        
        %% by diff max
        dif = diff(good_segments');
        dif = dif';
        good_seg_diff_peak = NaN(good_segment_count, 2);
        for I = 1:good_segment_count
            [good_seg_diff_peak(I, 1), good_seg_diff_peak(I, 2)] = max(dif(I, :));
        end


    else
        good_seg_diff_peak = NaN(good_segment_count, 2);
        for i = 1:size(good_segments,1)
            seg = good_segments(i,:);
%             half = round(length(seg)/2);
%             good_seg_diff_peak(:,1) = seg(half);
%             good_seg_diff_peak(:,2) = half;
            min_peak_height = nanmedian(seg)+ nanstd(seg) * 2.5;
            warning('off');
            [pks,locs] = findpeaks(seg,'MinPeakHeight',min_peak_height);
            warning('on');
            if isempty(locs) 
                good_seg_diff_peak(:,1) = seg(1);
                good_seg_diff_peak(:,2) = 1;
            else
                good_seg_diff_peak(:,1) = pks(1);
                good_seg_diff_peak(:,2) = locs(1);
            end
        end
    end
    
%     figure
%     s1 = good_segments(1,~isnan(good_segments(1,:)));
%     for i = 2:size(good_segments,1)
%         s2 = good_segments(i,~isnan(good_segments(i,:)));
%         X1=xcorr(s1,s2);   %compute cross-correlation between vectors s1 and s2
%         [~,d]=max(X1);      %find value and index of maximum value of cross-correlation amplitude
%         delay=d-max(length(s1),length(s2));   %shift index d, as length(X1)=2*N-1; where N is the length of the signals
%         hold on;
%         plot(delay+1:length(s2)+delay,s2,'r');   %Delay signal s2 by delay in order to align them
%     end
%     plot(s1);                             %Plot signal s1
%     grid on
%     
%     figure
%     s1 = good_segments(~isnan(good_segments(:,1)),1);
%     for i = 2:size(good_segments,1)
%         s2 = good_segments(~isnan(good_segments(:,i)),i);
%         X1=xcorr(s1,s2);   %compute cross-correlation between vectors s1 and s2
%         [~,d]=max(X1);      %find value and index of maximum value of cross-correlation amplitude
%         delay=d-max(length(s1),length(s2));   %shift index d, as length(X1)=2*N-1; where N is the length of the signals
%         hold on;
%         plot(delay+1:length(s2)+delay,s2,'r');   %Delay signal s2 by delay in order to align them
%     end
%     plot(s1);  
    alignedSegments = AlignWaves(good_segments, good_seg_diff_peak(:,2));
    alignment_point = max(good_seg_diff_peak(:, 2));
    
    %% mikes wave score
    compositeMean = NaN(1, size(alignedSegments, 2));
    compositeSTD = NaN(1, size(alignedSegments, 2));

    for I = 1:size(alignedSegments, 2)
        compositeMean(I) = nanmean(alignedSegments(:, I));
        compositeSTD(I) = nanstd(alignedSegments(:, I));
    end

    waveScore = zeros(1, size(alignedSegments, 1));

    for B = 1:size(alignedSegments, 1)

        for I = 1:size(alignedSegments, 2)

            if ((alignedSegments(B, I) > (compositeMean(I) + 2.5 * compositeSTD(I)) || alignedSegments(B, I) < (compositeMean(I) - 2.5 * compositeSTD(I))))
                waveScore(B) = waveScore(B) + 1;
            end

            if ((alignedSegments(B, I) > (compositeMean(I) + 3 * compositeSTD(I)) || alignedSegments(B, I) < (compositeMean(I) - 3 * compositeSTD(I))))
                waveScore(B) = waveScore(B) + 3;
            end

            if ((alignedSegments(B, I) > (compositeMean(I) + 4 * compositeSTD(I)) || alignedSegments(B, I) < (compositeMean(I) - 4 * compositeSTD(I))))
                waveScore(B) = waveScore(B) + 9;
            end

        end

    end

    waveScore = waveScore ./ good_seg_lengths';
    badamp = waveScore > .30;

    if dbg_plot == 1
        bad_type_count = [bad_type_count;numel(badamp(badamp == 1))];
        bad_total_count = bad_total_count + numel(badamp(badamp == 1));
        badSegments = alignedSegments(badamp == 1, :);
        UpdateBadSegFig(bad_segment_fig, badSegments, rgb('Salmon'));
    end
    
    good_segment_count = numel(badamp(badamp == 0));
    good_seg_lengths = good_seg_lengths(badamp == 0);
    good_segments = alignedSegments(badamp == 0, :);
    good_seg_diff_peak = good_seg_diff_peak(badamp == 0, :);


    %% Create Template
    beginpos = alignment_point - round(nanmedian(good_seg_diff_peak(:,2)));
    if beginpos < 1
        beginpos = 1;
    end
    endpos = alignment_point + round(nanmedian(good_seg_lengths - good_seg_diff_peak(:,2)));
    if endpos > size(good_segments,2)
        endpos = size(good_segments,2);
    end
    template = nanmedian(good_segments(:,beginpos:endpos));
    
    %% plot alignments
    if dbg_plot == 1
        
        UpdateBadSegFig(bad_segment_fig, good_segments, rgb('DodgerBlue'));

        subplot(3, 1, [1,2]);
        hold on;
        plot(nanmedian(good_segments)', 'LineWidth', 2, 'Color',rgb('Cyan'));
        
        MakeLegend(bad_segment_fig, {['Bad Length | ' num2str(bad_type_count(1))], ...
                                     ['Bad Peak amplitude | '  num2str(bad_type_count(2))], ...
                                     ['Bad Peak Position | ' num2str(bad_type_count(3))], ...
                                     ['Bad WaveScore | ' num2str(bad_type_count(4))], ...
                                     ['Good Waves | '  num2str(good_segment_count)], ...
                                     'Template'},...
                                    {rgb('black'),rgb('Lime'),rgb('Khaki'),rgb('Salmon'),rgb('DodgerBlue'),rgb('Cyan')});  
        
        vline(beginpos, {'m--', 'LineWidth', 2}, 'Template Begin');
        vline(endpos, {'m--','LineWidth', 2},'Template End');

        total_beats = size(segment_idxs,1);
        bad_precent = num2str((bad_total_count/total_beats) * 100, 3.5);
        good_precent = num2str((good_segment_count/total_beats) * 100, 3.5);
        title({'Press Space to continue.',['Total Beats: ' num2str(total_beats) ' | ' 'Bad: ' num2str(bad_total_count) ', ' bad_precent '% | '...
               'Good: ' num2str(good_segment_count) ', ' good_precent '%']});

        subplot(3, 1, 3);
        plot(template');
        title(fig_name);
        xlabel('Samples');
        ylabel('mV');
     end

end

function MakeLegend(fig,elements,colors)

    if length(elements) ~= length(colors)
       error('Element and color lengths must be same'); 
    end
    set(0, 'CurrentFigure', fig);
    subplot(3, 1, [1,2]);
    hold on;
    h = zeros(length(elements), 1);
    for i = 1:length(elements)
       h(i) = plot(NaN,NaN,'Color', colors{i});
    end
    xlabel('Samples');
    ylabel('mV');

    legend(h,elements);
    hold off;
end

function UpdateBadSegFig(axis, badSegs, color)
    set(0, 'CurrentFigure', axis);
    hold on;
    subplot(3, 1, [1,2]);
    p = plot(badSegs', 'Color', color);
    uistack(p,'bottom')
    hold off;
end