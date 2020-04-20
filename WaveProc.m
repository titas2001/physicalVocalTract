%{
    TODO: check virtual grid points, especially for boundary conditions
    NB: all equations except n°5 have wrong boundary conditions
%}

function [u,uNext] = WaveProc(uNext, u, uPrev, wNext, w, wPrev, lambdaSq, beta, k, h, N, L, c, S, rho, exciter, IR, number)
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
%exciter: single value, exciter signal
%IR: if !=0 bypasses exciter signal and excitates with impulse

% radiation parameters taken from Bilbao book
%alf = 2.0881*N*sqrt(1/(S(1)*S(N)));
% bet = 0.7407/c/N;
alf = L/(2*0.8216^2*c);
bet = L/(0.8216*sqrt(S(1)*S(2)/pi));

M = 0.01;
eps = c*sqrt(2*rho/M)*(pi/S(1))^(1/4);

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
         
     case 5  % Wave processing shape from Bilbao book no damp
         for l = 1:N
             if l == 1 %Closed end - Neumann condition: u(l-1) = u(l+1); S(l-1) = S(l+1)
                if (~IR)
                    shapeCoeff = (3*S(l) - S(l+1))/(2*S(l));
                    uNext(l) = 2*(1-lambdaSq)*u(l) - uPrev(l) + 2*lambdaSq * u(l+1) + (c^2*k^2/h)*exciter*shapeCoeff;
                else
                    uNext(l) = 2*lambdaSq*u(l+1) + 2*(1-lambdaSq)*u(l) - uPrev(l);
                end
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
        
     case 6  % Wave processing shape from Bilbao book damp
        for l = 1:N
            f0 = 80;
            Z = wPrev(l)*(beta*k-1) + w(l)*(1-k^2*f0^2) - k*eps*(S(l)^(1/4))*uPrev(l)/2;
            MAT1 = [- k*eps*(S(l)^(1/4))/2, 1+k*beta];
             if l == 1 %Closed end - Neumann condition: u(l-1) = u(l+1); S(l-1) = S(l+1)
                if (~IR)
                    shapeCoeff = (3*S(l) - S(l+1))/(2*S(l));
                    Q = 2*(1-lambdaSq)*u(l) - uPrev(l) + 2*lambdaSq * u(l+1) + (c^2*k^2/h)*exciter*shapeCoeff + eps*(S(l)^(1/4))*wPrev(l)/2;
                    MAT2 = [1, k*eps*(S(l)^(1/4))/2];
                else
                    Q = 2*lambdaSq*u(l+1) + 2*(1-lambdaSq)*u(l) - uPrev(l) + k*eps*(S(l)^(1/4))*wPrev(l)/(4*(S(l)+S(l+1)));
                    MAT2 = [1, k*eps*(S(l)^(1/4))/(4*(S(l)+S(l+1)))];
                end
             elseif l == N %Open end - loss condition
                coeff1 = (alf*h + bet*h*k)/k; %I radiating coefficient
                coeff2 = (alf*h - bet*h*k)/k; %II radiating coefficient
                coeff3 = 1+lambdaSq*coeff1;   %common denominator
                Q = lambdaSq*(2*u(l-1)+uPrev(l)*coeff2) + 2*(1-lambdaSq)*u(l) - uPrev(l) - eps*k*wPrev(l)*(S(l)^(1/4))/(4*S(l-1)+4*S(l));
                MAT2 = [coeff3, eps*k*(S(l)^(1/4))/(4*S(l-1)+4*S(l))];
             else   %Equation outside boundaries
                A = S(l+1) + 2*S(l) + S(l-1); %This arises from the double mean of S (Bilbao pg 256)
                coeff = (2*lambdaSq/A);       %Common coefficient
                Q = coeff*(S(l+1)+S(l))*u(l+1) + coeff*(S(l)+S(l-1))*u(l-1) + 2*(1-lambdaSq)*u(l) - uPrev(l) + k*eps*wPrev(l)*(S(l)^(1/4))/(2*A);
                MAT2 = [1, k*eps*(S(l)^(1/4))/(2*A)];
             end
             MAT = [MAT1; MAT2];
             MATinv = inv(MAT);
             paramVec = [Z; Q];
             resVec = MATinv*paramVec;
             uNext(l) = resVec(1);
             wNext(l) = resVec(2);
        end
end

end

