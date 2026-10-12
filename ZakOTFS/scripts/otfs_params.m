function p = otfs_params(M, N, B, D)
    p.M = M;  p.N = N;  p.B = B;  p.D = D;
    p.Ts        = 1/B;            % delay bin
    p.T         = M*p.Ts;         % one Zak period
    p.Tframe    = N*p.T;          % one OTFS frame
    p.Tsamp     = p.Ts/D;         % ADC/DAC sample interval
    p.fs        = 1/p.Tsamp;      % sample rate
    p.P         = M*D;            % samples per Zak period
    p.Lframe    = N*p.P;          % samples per frame
    p.frameRate = 1/p.Tframe;     % frames per second
    p.dtau      = p.Ts;           % delay resolution (slide p9)
    p.dnu       = 1/p.Tframe;     % Doppler spacing (derived from the DFT over N periods)
    p.nuSpan    = 1/p.T;          % total Doppler span (N*dnu)
end