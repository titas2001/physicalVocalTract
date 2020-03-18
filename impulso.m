function [imp] = impulso(mode, freq, fs, time)


% Inputs (change them at will)

% mode = refers to the type of input (thus far only a 
%               sine pulse and sawtooth)
% freq = frequency
% fs = sampling frequency
% time = length of the signal

t = 0:1/fs:time;

switch mode
    
    case 'sine_Pulse'
        
        imp = sin(2*pi*freq*t);
        imp(imp<0) = 0;
        
    case 'saw_Tooth'
        
        imp = sawtooth(2*pi*freq*t);
        imp(imp<0) = 0;

end