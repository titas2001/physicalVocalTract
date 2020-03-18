function [imp] = impulso(mode, freq, fs, time)


% Inputs (change them at will)

% mode = refers to the type of input (thus far only a 
%               sine pulse and sawtooth)
% freq = frequency
% fs = sampling frequency
% time = length of the signal

t = 1:time;

switch mode
    
    case 'sinepulse'
        
        imp = sin(2*pi*freq*t/fs);
        imp(imp<0) = 0;
        
    case 'sawpulse'
        
        imp = sawtooth(2*pi*freq*t/fs);
        imp(imp<0) = 0;
        
    case 'sawtooth'
        
        imp = sawtooth(2*pi*freq*t/fs);

end