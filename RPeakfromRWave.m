function ridxs = RPeakfromRWave(ecg,rWaveIdx)
    if length(rWaveIdx) <= 1
        ridxs = rWaveIdx;
        return;
    end
    avg_diff = median(diff(rWaveIdx));
    half_window_size = avg_diff/6;
%     close all;
% plot(rWaveIdx,ecg(rWaveIdx),'*r')
% hold on;
% plot(ecg)
    ridxs = rWaveIdx;
    for i = 1:length(rWaveIdx)
        if rWaveIdx(i)-half_window_size < 1
            winstart = 1; 
        else
            winstart = round(rWaveIdx(i)-half_window_size); 
        end
        if rWaveIdx(i)+half_window_size > length(ecg)
            winend = length(ecg); 
        else
            winend = round(rWaveIdx(i)+half_window_size); 
        end

        
        window = ecg(winstart:winend);
% vline(winstart)
% vline(winend)

        [~,id] = max(window);
        ridxs(i) = winstart+id-1;
%       plot(ridxs(i),ecg(ridxs(i)),'g*');
    end
% 
% 
%     for i = 1:length(rWaveIdx)
%         currpoint = rWaveIdx(i);
%         while 1
%             if currpoint > 1 && currpoint < length(ecg)
%                 if ecg(currpoint) <= ecg(currpoint+1)
%                     currpoint = currpoint +1;
%                 elseif ecg(currpoint) < ecg(currpoint-1)
%                     currpoint = currpoint - 1;
%                 else
%                     ridxs(i) = currpoint;
%                     break
%                 end
%             else
%                 ridxs(i) = currpoint;
%                 break
%             end
% 
%         end
%     end
%     time = 1:length(ecg);
%     plot(time,ecg);
%     hold on;
%     scatter(time(rWaveIdx),ecg(rWaveIdx),'r');
%     scatter(time(ridxs),ecg(ridxs),'c');


end