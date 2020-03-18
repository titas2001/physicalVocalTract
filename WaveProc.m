function [u,uNext] = WaveProc(uNext, u, uPrev, lambdaSq, beta, k, N, S, number)
% Wave Processing: runs selected update equation call function:
% [u,uNext] = WaveProc(uNext, u, uPrev, lambdaSq, beta, k, N, S, number)
% choose equation with 'number':
%    0 --> Wave processing damp backwards derivative
%    1 --> Wave processing damp center derivative
%    2 --> Wave processing shape damp center derivative
%    3 --> Wave processing no damp shape function
%    4 --> Wave processing no damp

switch number
    case 0
        % Wave processing damp backwards derivative
         for l = 2:N
             if l == N %Free end
                uNext(l) = (2-beta*k) * u(l) + (beta*k-1) * uPrev(l) + lambdaSq * (2 * u(l-1) - 2 * u(l)); 
             else
                uNext(l) = (2-beta*k) * u(l) + (beta*k-1) * uPrev(l) + lambdaSq * (u(l+1) - 2 * u(l) + u(l-1));
             end
         end
    case 1    
        % Wave processing damp center derivative
         for l = 2:N
             if l == N %Free end
                uNext(l) = (4/(2+beta*k)) * u(l) + ((beta*k-2)/(beta*k+2)) * uPrev(l) + (2*lambdaSq/(2+beta*k)) * (2 * u(l-1) - 2 * u(l)); 
             else
                uNext(l) = (4/(2+beta*k)) * u(l) + ((beta*k-2)/(beta*k+2)) * uPrev(l) + (2*lambdaSq/(2+beta*k)) * (u(l+1) - 2 * u(l) + u(l-1));
             end
         end
    case 2
        % Wave processing shape damp center derivative
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
    case 3
        % Wave processing no damp shape function
         for l = 2:N
             Smean = (S(l) + S(l+1))/2;
             if l == N %Free end
                uNext(l) = 2*(1-lambdaSq) * u(l) - uPrev(l) + (lambdaSq * S(l+1)/Smean) * u(l-1) + (lambdaSq * S(l)/Smean) * u(l-1); 
             else
                uNext(l) = 2*(1-lambdaSq) * u(l) - uPrev(l) + (lambdaSq * S(l+1)/Smean) * u(l+1) + (lambdaSq * S(l)/Smean) * u(l-1);
             end
         end
    case 4
        % Wave processing no damp
         for l = 2:N
             if l == N %Free end
                uNext(l) = 2 * u(l) - uPrev(l) + lambdaSq * (2 * u(l-1) - 2 * u(l)); 
             else
                uNext(l) = 2 * u(l) - uPrev(l) + lambdaSq * (u(l+1) - 2 * u(l) + u(l-1));
             end
         end
end

end

