function pulseQualityIndex = PQI(rawPPGsignal, samplingRate)
    
    % Call Pre-Processing
    [PPG25dHz, PPG10Hz] = preprocessing(rawPPGsignal);
    
    % call (...)
    [PpTemp, amplitudeCorrectionFactors, PpPQI] = PPGsegmentationAndBeatLocalization(PPG25dHz, PPG10Hz);
    
    % call (...)
    % per each pulse -> loop?
    templateAD = templateCreartion(PpTemp, amplitudeCorrectionFactors, PpPQI);
    
    % call (...)
    pulseQualityIndex = pulseTemplateComparision(templateAD);
    
end
