function peaks = hill_climber(epoch,pan_tomp_Rpeaks,sample_rate)


% included 45, as wide qrs in 0.12s, take three times that, 0.36s and
% convert to samples
% but in future do hill climbing
peaks = [];
search_Space = ceil((25*sample_rate)/125);
for i = 1:length(pan_tomp_Rpeaks)
    i;
    current_R_peak = pan_tomp_Rpeaks(i);
    search_Space_after_peak = ceil(current_R_peak + search_Space);
    search_Space_before_peak = floor(current_R_peak - search_Space);
    if search_Space_after_peak>length(epoch)
        search_Space_after_peak = length(epoch);
       [~,max_peak] = max(epoch(search_Space_before_peak:search_Space_after_peak));
       peaks = [peaks,max_peak+search_Space_before_peak];

    elseif search_Space_before_peak<1
        search_Space_before_peak = 1;
       [~,max_peak] = max(epoch(search_Space_before_peak:search_Space_after_peak));
       peaks = [peaks,max_peak+search_Space_before_peak];

    else
        [~,max_peak] = max(epoch(search_Space_before_peak:search_Space_after_peak));
        max_peak = search_Space_before_peak + max_peak -1;
        peaks = [peaks,max_peak];
    end
    
end