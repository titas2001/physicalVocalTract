%{
    TODO: check virtual grid points, especially for boundary conditions
    NB: all equations except n°5 have wrong boundary conditions
%}

function [u,uNext] = WaveProc(uNext, u, uPrev, lambdaSq, beta, k, h, N, L, c, S, number, bound)
% Wave Processing: runs selected update equation call function:
% [u,uNext] = WaveProc(uNext, u, uPrev, lambdaSq, beta, k, N, S, number)
% choose equation with 'number':
%    0 --> Wave processing damp backwards derivative
%    1 --> Wave processing damp center derivative
%    2 --> Wave processing shape damp center derivative
%    3 --> Wave processing no damp shape function
%    4 --> Wave processing no damp
%    5 --> Wave processing shape from Bilbao book
%
% bound: bounding condition at open end (left)
%    0 --> Dirichlet (not working)
%    1 --> Energy loss condition

% radiation parameters taken from Bilbao book
%alf = 2.0881*N*sqrt(1/(S(1)*S(N)));
% bet = 0.7407/c/N;
alf = L/(2*0.8216^2*c);
bet = L/(0.8216*sqrt(S(1)*S(2)/pi));

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
         
     case 5  % Wave processing shape from Bilbao book
        for l = 1:N
             if l == 1 %Closed end - Neumann condition: u(l-1) = u(l+1); S(l-1) = S(l+1)
                uNext(l) = 2*lambdaSq*u(l+1) + 2*(1-lambdaSq)*u(l) - uPrev(l);
             elseif l == N %Open end - loss condition
                coeff1 = (alf*h + bet*h*k)/k; %I radiating coefficient
                coeff2 = (alf*h - bet*h*k)/k; %II radiating coefficient
                coeff3 = 1+lambdaSq*coeff1;   %common denominator
                uNext(l) = (lambdaSq/coeff3)*(2*u(l-1)+uPrev(l)*coeff2) + (2*(1-lambdaSq)*u(l) - uPrev(l))/coeff3;
             else   %Equation outside boundaries
                A = S(l+1) + 2*S(l) + S(l-1); %This arises from the double mean of S (Bilbao pg 256)
                coeff = (2*lambdaSq/A);       %Common coefficient
                uNext(l) = coeff*(S(l+1)+S(l))*u(l+1) + coeff*(S(l)+S(l-1))*u(l-1) + 2*(1-lambdaSq)*u(l) - uPrev(l);
             end
        end
end

end

