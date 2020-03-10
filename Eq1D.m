clear all; 
close all;

fs = 44100;         % sample rate
k = 1 / fs;         % time step
dur = fs;           % duration

% Define variables of the system
c = 344;

% Calculate grid spacing from variables
h = c * k;
N = floor(1/h);
%h = 3 / N;

% Calculate courant number
lambdaSq = c^2 * k^2 / h^2

% Initialise u^{n+1} and u^n
uNext = zeros(N, 1);
u = zeros(N, 1);
u(2) = 1; %exiciting the first samples
u(1) = 0; 
% initialise u^{n-1}
uPrev = u;

% Initialise output and output position
out = zeros(dur, 1);
outPos = N;

%% Tube Function

S = ones(N,1);
D = N - floor(N/4);
for i = D:N
    S(i) = exp((i-D)*0.04);
end
%plot(S)

exciter = zeros(1,dur);

for n = 1:dur 
     %Exciter
     exciter(n) = sawtooth(2*pi*20*(n-1)/fs);
     %u(1) = exciter(n);
     
     for l = 2:N
         if l == N %Free end
            uNext(l) = 2 * u(l) - uPrev(l) + lambdaSq * (2 * u(l-1) - 2 * u(l)); 
         else
            uNext(l) = 2 * u(l) - uPrev(l) + lambdaSq * (u(l+1) - 2 * u(l) + u(l-1));
         end
     end
    % retrieve output
    out(n) = uNext(outPos);
    
    % draw string
    %plot(uNext);
    %plot(exiciter);
    %ylim([-1, 1]);
    %drawnow;
    
    uPrev = u;
    u = uNext;
end

soundsc(out, fs);

%Plotting Output
freqScaling = fs/dur;
freqAxis = freqScaling:freqScaling:(freqScaling*dur);
transform = abs(fft(out));
tiledlayout(2,1)
% Top plot
nexttile
plot(out)
title('Time')
% Bottom plot
nexttile
plot(freqAxis(1:22050),transform(1:22050))
title('Freq')

%plot(exiciter)