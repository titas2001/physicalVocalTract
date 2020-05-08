clear all; 
close all;

%If looking for Impulse Response set IR != 0
IR = 0;

%Shape parameters
curveStartPos = 0.8;
maxWidth = 0.06;
minWidth = 0.01;
shapeType = 'aeB';

%fs = 44100;                % time grid sample rate
%fs = 157870;               % time grid sample rate
fs = 480e3;                 % time grid sample rate
Fs = 48e3;                  % playback sample rate
k = 1 / fs;                 % time step [s]
durSec = 3;                 % duration in seconds
dur = durSec*fs;                 % duration [samples]
playbackDur = Fs*durSec;

%{ 
    scaling factor used to downsample before playback and plotting
    it has to be an integer!
%}
sFactor = fs/Fs;            


%Damp coefficient
beta = 8125;

% Define speed of sound - Assume 37°C c~331.6+0.6*T
c = 353.8;                    %[m/s]

% Air density at 37°C, 1 atm
rho = 1.138; 

%Mass per unit length given epsilon value from Bilbao
M = 0.01; %0.3947;
%M=4;
M = 0.012; %0.3947;
%resonant frequency, from Bilbao
f0 = 0;

% Calculate grid spacing from variables
h = c * k;
% L = 0.1765;                      % Tube length [m] 
% L = 0.165475;

if shapeType(1) == 'u'
    L = 0.195;
else
    L = 0.165475;
end

% L = 1;
N = floor(L/h);             % Tube length [samples]
h = L / N;
d = (c*sqrt(N))/fs;
% Calculate courant number
lambdaSq = c^2 * k^2 / h^2;

% Exciter frequency
exciterFreq = 120;           % [Hz]

% Initialise spatial states u(n+1) and u(n)
uNext = zeros(N, 1);
u = zeros(N, 1);

wNext = zeros(N, 1);
w = zeros(N, 1);
wPrev = zeros(N, 1);

% Exiciting with impulse at closed end
% (Impulse response)
if IR
    %u(2) = 1;              % raw impulse
    width = floor(N/10);    % hann impulse 
    u(1:width) = hann(width);
end

% Initialise spatial state u(n-1)
uPrev = u;

% Initialise output vector
out = zeros(playbackDur, 1);
outUpSample = zeros(fs*durSec,1);

% Defining where output is observed, in our case the end of the tube
outPos = N;

%Shape function (*2 because we look for the area, not sure about this)
S = Shape(N+1, curveStartPos, minWidth, maxWidth, shapeType)*2;

exciterSign = Impulso('camel', exciterFreq, fs, dur, 45);

fil = VT_Filter(fs, 16000, 20000);
numnum = fil.Numerator;

for n = 1:dur
    
    %TODO: Bilbao does this, have to understand why
    excit = 0.5*(exciterSign(n)+abs(exciterSign(n)));
    %excit = exciterSign(n);
    [u,uNext] = WaveProc(uNext, u, uPrev, wNext, w, wPrev, lambdaSq, beta, k, h, N, L, c, S, rho, M, f0, excit, IR, 6);
    
    %     Retrieve output, p=(c^2ro/S)dphi/dt, filling output vector
    %outUpSample(n) = (rho*c^2/S(N))*(uNext(outPos) - u(outPos)) / k;
    outUpSample(n) = rho*(uNext(outPos) - u(outPos)) / k;
    
% %     Retrieve output, p=(c^2ro/S)dphi/dt, filling output vector
%     if mod(n, sFactor) == 0
%         out(counter) = VT_fil(n);
%         counter = counter + 1;
%     end
%     %out(n) = uNext(outPos);

%     % Real time states drawing
%     plot(uNext);
%     %plot(exiciter);
%     ylim([-0.01, 0.01]);
%     drawnow;
    
    % Update spatial states
    uPrev = u;
    u = uNext;
end

VT_fil = filter(numnum, 1, outUpSample);
counter = 1;
for i=1:size(outUpSample)
    %     Retrieve output, p=(c^2ro/S)dphi/dt, filling output vector
    if mod(i, sFactor) == 0
        out(counter) = VT_fil(i);
        counter = counter + 1;
    end
end

% Normalizing output
maxOut = max(out);    % find max value of output
minOut = abs(min(out));
if maxOut>minOut
    out = out/maxOut;
else
    out = out/minOut;
end
% i = 1:floor(length(out)/sFactor);
% nOut = out(i)/mVal;    % normalized Output

%nOut = lowpass(out, 0.0119*sFactor);
nOut = out;
sound(nOut, Fs);

%audiowrite("French_U_16_5.wav",nOut,Fs);

%Plotting Output
freqScaling = Fs/playbackDur;
freqAxis = freqScaling:freqScaling:(freqScaling*playbackDur);
transform = abs(fft(nOut));
figure(1)
tiledlayout(3,1)
% Top plot
nexttile
plot(nOut)
title('Time')
xlabel('samples')
% Bottom plot
nexttile
plot(freqAxis(1:(66150)),transform(1:(66150)))
title('Freq')
xlabel('Hz')
nexttile
plot(S/2);
hold on
plot(-S/2);
title('Shape')

figure(2)
plot(exciterSign)