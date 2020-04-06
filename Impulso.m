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
        
    case 'camel'
        
        t = 0.001:0.001:1;
        e1 = 0.1;
        e2 = 0.4;
        
        for i=1:length(t)
            if t(i)<e1
               x(i) = 0.5-0.5*cos(2*pi*freq*t(i));
            elseif t(i)<=e2
               x(i) = -((0.5-0.5*cos(2*pi*freq*e1))/(e2-e1))*((t(i)-e2));
            else
               x(i) = 0;
            end
        end
         xx = repmat(x, [1, time]);
         imp = xx;

end