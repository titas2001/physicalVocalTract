clear all; 
close all;

%If looking for Impulse Response set IR != 0
IR = 0;

%Shape parameters
curveStartPos = 0.8;
maxWidth = 0.06;
minWidth = 0.01;
shapeType = 'exp';

%Damp coefficient
beta = 0.3;

fs = 44100;         % sample rate
k = 1 / fs;         % time step
dur = 3*fs;           % duration

% Define speed of sound
c = 344;

% Calculate grid spacing from variables
h = c * k;
N = floor(1/h); %length of the tube
%h = 1 / N;

% Calculate courant number
lambdaSq = c^2 * k^2 / h^2;

% Initialising excitation function
exciter = zeros(1,dur);
exticerFreq = 100;

% Initialise spatial states u(n+1) and u(n)
uNext = zeros(N, 1);
u = zeros(N, 1);

% Exiciting with impulse at closed end
% (Impulse response)
if IR
    %u(2) = 1; 
    width = floor(N/10);
    u(1:width) = hann(width); %More physical impulse
end

% Initialise spatial state u(n-1)
uPrev = u;

% Initialise output vector
out = zeros(dur, 1);

% Defining where output is observed, in our case the end of the tube
outPos = N;

%Shape function (*2 because we loor for the area, not sure about this
%though)
S = Shape(N+1, curveStartPos, minWidth, maxWidth, shapeType) * 2;

for n = 1:dur 
     %Exciter processing
     exciter(n) = sawtooth(2*pi*exticerFreq*(n-1)/fs);
     if ~IR
         %Multiplying by shape so it's not [-1,1] because it seems
         %reasonable but maybe it's not
        u(2) = exciter(n) * S(2)/2;
     end
     [u,uNext] = WaveProc(uNext, u, uPrev, lambdaSq, beta, k, N, S, 2);
     
    % Retrieve output, filling output vector
    out(n) = uNext(outPos);
    
%     % Real time states drawing
%     plot(uNext);
%     %plot(exiciter);
%     ylim([-1, 1]);
%     drawnow;
    
    % Update spatial states
    uPrev = u;
    u = uNext;
end

soundsc(out, fs);

%Plotting Output
freqScaling = fs/dur;
freqAxis = freqScaling:freqScaling:(freqScaling*dur);
transform = abs(fft(out));
figure(1)
tiledlayout(3,1)
% Top plot
nexttile
plot(out)
title('Time')
% Bottom plot
nexttile
plot(freqAxis(1:66150),transform(1:66150))
title('Freq')
nexttile
plot(S);
title('Shape')

%plot(exiciter)