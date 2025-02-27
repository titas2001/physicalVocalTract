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
pressure = 0.09548; %~70 dB
switch mode
    
    case 'sinepulse'
        
        imp = sin(2*pi*freq*t/fs) - 1 + (dutycycle/100)*2;
        imp(imp<0) = 0;
        imp = imp/max(imp);
        imp=imp*pressure;
        
    case 'sawpulse'
        
        imp = sawtooth(2*pi*freq*t/fs) - 1 + (dutycycle/100);
        imp(imp<0) = 0;
        imp = imp/max(imp);
        imp=imp*pressure;
        
    case 'sawtooth'
        
        imp = sawtooth(2*pi*freq*t/fs);
        imp=imp*pressure;
        
    case 'glottal_pulse'
        
        t = freq/fs:freq/fs:1+freq/fs;
        e1 = 0.1;
        e2 = 0.4;
        repeat = 3*freq;
        
        for i=1:length(t)
            if t(i)<e1
               x(i) = 0.5-0.5*cos(2*pi*t(i));
            elseif t(i)<=e2
               x(i) = -((0.5-0.5*cos(2*pi*e1))/(e2-e1))*((t(i)-e2));
            else
               x(i) = 0;
            end
        end

         xx = repmat(x, [1, repeat]);
         wit_noise = awgn(xx, -24, 'measured');
         oioi = xx+0.0001*wit_noise;
         imp = oioi;
         imp = imp/max(imp);
         imp=imp*pressure;
end