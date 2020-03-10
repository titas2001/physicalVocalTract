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
h = 1 / N;

% Calculate courant number
lambdaSq = c^2 * k^2 / h^2;

% Initialise u^{n+1} and u^n
uNext = zeros(N, 1);
u = zeros(N, 1);
%u(2) = 1; %exiciting the first sample
% initialise u^{n-1}
uPrev = u;

% Initialise output and output position
out = zeros(dur, 1);
outPos = N-1;

%% Loop
%range = 2:N-1; % define range for which we want to calculate (boundary condition)

exciter = zeros(1,dur);

for n = 1:dur 
     %nested for-loop
     exciter(n) = sawtooth(2*pi*100*n/fs);
     u(1) = exciter(n);
     %u(1)
     for l = 2:N
         if l == N %Free end
            uNext(l) = 2 * u(l) - uPrev(l) + lambdaSq * (2 * u(l-1) - 2 * u(l)); 
         else
            uNext(l) = 2 * u(l) - uPrev(l) + lambdaSq * (u(l+1) - 2 * u(l) + u(l-1));
         end
         %u(1)
     end
    
    % retrieve output
    out(n) = uNext(outPos);
    
    % draw string
    plot(uNext);
    %plot(exiciter);
    ylim([-1, 1]);
    drawnow;
    
    uPrev = u;
    u = uNext;
end

%sound(out, fs);
%plot(out)
%plot(exiciter)