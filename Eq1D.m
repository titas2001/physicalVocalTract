clear all; 
close all;

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
lambdaSq = c^2 * k^2 / h^2

%Damp coefficient
%(something weird happens if beta>=3.427, investigate)
beta = 2;

% Initialise spatial states u(n+1) and u(n)
uNext = zeros(N, 1);
u = zeros(N, 1);

% Exiciting with impulse at closed end
% (Impulse response)
%u(2) = 1; 

% Initialise spatial state u(n-1)
uPrev = u;

% Initialise output vector
out = zeros(dur, 1);

% Defining where output is observed, in our case the end of the tube
outPos = N;

% Initialising excitation function
exciter = zeros(1,dur);

for n = 1:dur 
     %Exciter processing
     exciter(n) = sawtooth(2*pi*80*(n-1)/fs);
     u(2) = exciter(n);
     
     %Wave processing damp backwards derivative
     for l = 2:N
         if l == N %Free end
            uNext(l) = (2-beta*k) * u(l) + (beta*k-1) * uPrev(l) + lambdaSq * (2 * u(l-1) - 2 * u(l)); 
         else
            uNext(l) = (2-beta*k) * u(l) + (beta*k-1) * uPrev(l) + lambdaSq * (u(l+1) - 2 * u(l) + u(l-1));
         end
     end
     
%      %Wave processing damp center derivative
%      for l = 2:N
%          if l == N %Free end
%             uNext(l) = (4/(2+beta*k)) * u(l) + (beta*beta*k*k/4 - 1) * uPrev(l) + (2*lambdaSq/(2+beta*k)) * (2 * u(l-1) - 2 * u(l)); 
%          else
%             uNext(l) = (4/(2+beta*k)) * u(l) + (beta*beta*k*k/4 - 1) * uPrev(l) + (2*lambdaSq/(2+beta*k)) * (u(l+1) - 2 * u(l) + u(l-1));
%          end
%      end

%      %Wave processing no damp
%      for l = 2:N
%          if l == N %Free end
%             uNext(l) = 2 * u(l) - uPrev(l) + lambdaSq * (2 * u(l-1) - 2 * u(l)); 
%          else
%             uNext(l) = 2 * u(l) - uPrev(l) + lambdaSq * (u(l+1) - 2 * u(l) + u(l-1));
%          end
%      end
     
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