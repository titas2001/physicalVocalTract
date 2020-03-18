clear all; 
close all;

%If looking for Impulse Response set IR != 0
IR = 0;

%Shape parameters
curveStartPos = 0.6;
maxWidth = 1;
minWidth = 0.02;
shapeType = 'linear';

%Damp coefficient
beta = 1;

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
        u(2) = exciter(n);
     end
     
%      %Wave processing damp backwards derivative
%      for l = 2:N
%          if l == N %Free end
%             uNext(l) = (2-beta*k) * u(l) + (beta*k-1) * uPrev(l) + lambdaSq * (2 * u(l-1) - 2 * u(l)); 
%          else
%             uNext(l) = (2-beta*k) * u(l) + (beta*k-1) * uPrev(l) + lambdaSq * (u(l+1) - 2 * u(l) + u(l-1));
%          end
%      end
     
% %      Wave processing damp center derivative
%      for l = 2:N
%          if l == N %Free end
%             uNext(l) = (4/(2+beta*k)) * u(l) + ((beta*k-2)/(beta*k+2)) * uPrev(l) + (2*lambdaSq/(2+beta*k)) * (2 * u(l-1) - 2 * u(l)); 
%          else
%             uNext(l) = (4/(2+beta*k)) * u(l) + ((beta*k-2)/(beta*k+2)) * uPrev(l) + (2*lambdaSq/(2+beta*k)) * (u(l+1) - 2 * u(l) + u(l-1));
%          end
%      end

%Wave processing shape damp center derivative
     for l = 2:N
         Smean = (S(l) + S(l+1))/2;
         coeff = 2*Smean+beta*k;
         if l == N %Free end
            spacePart = S(l+1)*(u(l) - u(l-1)) + S(l)*(3*u(l-1) - 3*u(l));
         else
            spacePart = S(l+1)*(u(l) - u(l-1)) + S(l)*(u(l+1) - 3*u(l) + 2*u(l-1));
         end
         uNext(l) = (4*Smean/coeff) * u(l) + ((beta*k-2*Smean)/coeff) * uPrev(l) + (2*lambdaSq/coeff) * spacePart;
     end

% %      Wave processing no damp shape function
%      for l = 2:N
%          Smean = (S(l) + S(l+1))/2;
%          if l == N %Free end
%             uNext(l) = 2*(1-lambdaSq) * u(l) - uPrev(l) + (lambdaSq * S(l+1)/Smean) * u(l-1) + (lambdaSq * S(l)/Smean) * u(l-1); 
%          else
%             uNext(l) = 2*(1-lambdaSq) * u(l) - uPrev(l) + (lambdaSq * S(l+1)/Smean) * u(l+1) + (lambdaSq * S(l)/Smean) * u(l-1);
%          end
%      end

% %      Wave processing no damp
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