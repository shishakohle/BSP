%% Function PQI() as a solution to Assignment 3

% expected lag between R peak in ECG to PPG beat is the pule transit time (PTT) fro the heart to the tip of teh finger ~200ms SOURCE?

function real_beats = findbeatsfromECGinPPG(ECGtimes, PPGtimes)
    
    PTT = 0.2; % time in ms
    ECGplusPTTtimes = ECGtimes + PTT;

%     for i = 1 : size(ECGtimes, 1)
        
        real_beats = find(PPGtimes >= (ECGplusPTTtimes - 0.05) | PPGtimes <= (ECGplusPTTtimes + 0.05));
        
%     end

end

% validate beat detection - ToDo Ingo

% validate inter beat interval durations - ToDo Alex