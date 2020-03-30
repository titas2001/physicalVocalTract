%{
TODO: check virtual grid points, especially for boundary conditions
%}

function [u,uNext] = WaveProc(uNext, u, uPrev, lambdaSq, beta, k, h, N, c, S, number, bound)
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

% radiation parameters
alf = 2.0881*N*sqrt(1/(S(1)*S(N)));
bet = 0.7407/c/N;

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
     case 5
    % Wave processing shape from Bilbao book
        for l = 1:N
             if l == 1 %Closed end: u(l-1) = u(l); S(l-1) = S(l)
                A = S(l+1) + 3*S(l);
                coeff = (2*lambdaSq/A);
                uNext(l) = coeff*(S(l+1)+S(l))*u(l+1) + 2*(coeff*S(l) + 1-lambdaSq)*u(l) - uPrev(l); 
             elseif l == N %Open end
                    A = 3*S(l) + S(l-1);
                    coeff = (2*lambdaSq/A);
                    if bound == 0 %dirichlet condition
                        uNext(l) = (2*coeff*S(l) + 2*(coeff*S(l) + 1-lambdaSq))*u(l) - uPrev(l); 
                    elseif bound == 1 %loss condition
                        %backwards derivative
        %                     coeff1 = 2*k + 2*alf*h + k*bet*h;
        %                     coeff2 = 2*alf*h - k*bet*h;
        %                     boundTerm = (2*k/coeff1) * u(l-1) + (coeff2/coeff1)*uPrev(l);
        %                     uNext(l) = (coeff*(S(l+1)+S(l)) + 2*(coeff*S(l) + 1-lambdaSq))*boundTerm - uPrev(l);
                        %center derivative
                        coeff1 = (alf*h + bet*h*k)/(2*k);
                        coeff2 = (alf*h - bet*h*k)/(2*k);
                        coeff3 = 1+2*S(l)*coeff*coeff1;
                        uNext(l) = (2*coeff*S(l)/coeff3)*(u(l)+uPrev(l)*coeff2) + (coeff*(S(l)+S(l-1))*u(l-1) + 2*(1-lambdaSq)*u(l) - uPrev(l))/coeff3;
                    end
             else
                A = S(l+1) + 2*S(l) + S(l-1);
                coeff = (2*lambdaSq/A);
                uNext(l) = coeff*(S(l+1)+S(l))*u(l+1) + coeff*(S(l)+S(l-1))*u(l-1) + 2*(1-lambdaSq)*u(l) - uPrev(l);
             end
        end
end

end

