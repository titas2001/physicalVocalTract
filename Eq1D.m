clear all; 
close all;

%If looking for Impulse Response set IR != 0
IR = 0;

%Shape parameters
curveStartPos = 0.8;
maxWidth = 0.06;
minWidth = 0.01;
shapeType = 'a';

%Damp coefficient (not used yet)
beta = 0.3;

% fs = 44100;                 % sample rate
fs = 157870;                 % sample rate
k = 1 / fs;                 % time step [s]
dur = 3*fs;                 % duration [samples]

% Define speed of sound
c = 344;                    %[m/s]

% Air density at 15°C, 1 atm
rho = 1.115; 

% Calculate grid spacing from variables
h = c * k;
L = 0.1765;                      % Tube length [m] 
% L = 1;
N = floor(L/h);             % Tube length [samples]
h = L / N;
d = (c*sqrt(N))/fs;
% Calculate courant number
lambdaSq = c^2 * k^2 / h^2;

% Exciter frequency
exticerFreq = 50;           % [Hz]

% Initialise spatial states u(n+1) and u(n)
uNext = zeros(N, 1);
u = zeros(N, 1);

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
out = zeros(dur, 1);

% Defining where output is observed, in our case the end of the tube
outPos = N;

%Shape function (*2 because we look for the area, not sure about this)
S = Shape(N+1, curveStartPos, minWidth, maxWidth, shapeType)*2;

%Initializing exciter
exciter = Impulso('sinepulse', exticerFreq, fs, dur, 45);

for n = 1:dur
    if ~IR
        %TODO: change the way u is excited (Bilbao book)
       u(2) = exciter(n) * S(2)/2;
    end
    [u,uNext] = WaveProc(uNext, u, uPrev, lambdaSq, beta, k, h, N, L, c, S, 5, 1);
     
    % Retrieve output, p=(c^2ro/S)dphi/dt, filling output vector
    out(n) = (rho*c^2/S(N))*(uNext(outPos) - u(outPos)) / k;
    %out(n) = uNext(outPos);

    %     % Real time states drawing
%         plot(uNext);
%         %plot(exiciter);
%         ylim([-0.01, 0.01]);
%         drawnow;
    
    % Update spatial states
    uPrev = u;
    u = uNext;
end
% y = downsample(out,10);
% soundsc(out, fs);

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
xlabel('samples')
% Bottom plot
nexttile
plot(freqAxis(1:66150),transform(1:66150))
title('Freq')
xlabel('Hz')
nexttile
plot(S/2);
hold on
plot(-S/2);
title('Shape')

% figure(2)
% plot(exciter)