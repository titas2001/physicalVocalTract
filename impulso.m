function [imp] = Impulso(mode, freq, fs, time, dutycycle)

% Parameters

% [mode]        determines type of signal output:
%                   'sinepulse' pulse train made from sinewave
%                   'sawpulse'  pulse train made from sawtooth
%                   'sawtooth'  simple test sawtooth wave
% [freq]        frequency of the signal (int>0)
% [fs]          sampling frequency (int>0)
% [time]        length of the signal (int>0)
% [dutycycle]   [1,100] in case of a pulse train determines fraction of
%               period each pulse occupies, i.e. 50 is a regular clipped
%               wave

if dutycycle<1 || dutycycle>100
        error('ducycycle has to be an int between 1 and 100');
end

if freq < 0
        error('freq has to be an int > 0');
end

if fs < 0
        error('fs has to be an int > 0');
end

if time < 0
        error('time has to be a scalar int > 0');
end

t = 1:time;

switch mode
    
    case 'sinepulse'
        
        imp = sin(2*pi*freq*t/fs) - 1 + (dutycycle/100)*2;
        imp(imp<0) = 0;
        imp = imp/max(imp);
        
    case 'sawpulse'
        
        imp = sawtooth(2*pi*freq*t/fs) - 1 + (dutycycle/100);
        imp(imp<0) = 0;
        imp = imp/max(imp);
        
    case 'sawtooth'
        
        imp = sawtooth(2*pi*freq*t/fs);

end