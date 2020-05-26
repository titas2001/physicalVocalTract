function [u,uNext] = WaveProc(uNext, u, uPrev, wNext, w, wPrev, lambdaSq, beta, k, h, N, L, c, S, rho, M, f0, exciter, IR, number)
% Wave Processing: runs selected update equation call function:
% choose equation with 'number' (number 2 is deprecated)
%exciter: single value, exciter signal
%IR: if !=0 bypasses exciter signal and excitates with impulse

% radiation parameters taken from Bilbao book
alf = 1/(2*0.8216^2*c);
bet = L/(0.8216*sqrt(S(1)*S(N)/pi));

eps = c*sqrt(2*rho/M)*(pi/S(1))^(1/4);

switch number
    case 1  % Wave processing shape from Bilbao book damp
        for l = 1:N
            Z = wPrev(l)*(beta*k-1) + w(l)*(1-k^2*f0^2) - k*eps*(S(l)^(1/4))*uPrev(l)/2;
            mat11 = - k*eps*(S(l)^(1/4))/2;
            mat12 =  1+k*beta;
             if l == 1 %Closed end - Neumann condition: u(l-1) = u(l+1); S(l-1) = S(l+1)
                if (~IR)
                    shapeCoeff = (3*S(l) - S(l+1))/(2*S(l));
                    Q = 2*(1-lambdaSq)*u(l) - uPrev(l) + 2*lambdaSq * u(l+1) + (c^2*k^2/h)*exciter*shapeCoeff + 2*k*eps*(S(l)^(1/4))*wPrev(l)/S(l);
                    mat21 = 1;
                    mat22 = 2*k*eps*(S(l)^(1/4))/S(l);
                else
                    Q = 2*lambdaSq*u(l+1) + 2*(1-lambdaSq)*u(l) - uPrev(l) + k*eps*(S(l)^(1/4))*wPrev(l)/(4*(S(l)+S(l+1)));
                    mat21 = 1;
                    mat22 = k*eps*(S(l)^(1/4))/(4*(S(l)+S(l+1)));
                end
             elseif l == N %Open end - loss condition
                coeff1 = (alf*h + bet*h*k)/k; %I radiating coefficient
                coeff2 = (alf*h - bet*h*k)/k; %II radiating coefficient
                coeff3 = 1+lambdaSq*coeff1;   %common denominator
                Q = lambdaSq*(2*u(l-1)+uPrev(l)*coeff2) + 2*(1-lambdaSq)*u(l) - uPrev(l) - eps*k*wPrev(l)*(S(l)^(1/4))/(4*S(l-1)+4*S(l));
                mat21 = coeff3;
                mat22 = eps*k*(S(l)^(1/4))/(4*S(l-1)+4*S(l));
             else   %Equation outside boundaries
                A = S(l+1) + 2*S(l) + S(l-1); %This arises from the double mean of S (Bilbao pg 256)
                coeff = (2*lambdaSq/A);       %Common coefficient
                Q = coeff*(S(l+1)+S(l))*u(l+1) + coeff*(S(l)+S(l-1))*u(l-1) + 2*(1-lambdaSq)*u(l) - uPrev(l) + k*eps*wPrev(l)*(S(l)^(1/4))/(2*A);
                mat21 = 1;
                mat22 = k*eps*(S(l)^(1/4))/(2*A);
             end
             %Cramer rule for 2 equations systems
             detMat = mat11*mat22 - mat12*mat21;
             detMatPsi = Z*mat22 - mat12*Q;
             detMatW = mat11*Q - Z*mat21;

             uNext(l) = detMatPsi/detMat;
             wNext(l) = detMatW/detMat;
        end
        
     case 2  % DEPRECATED Wave processing shape from Bilbao book no damp 
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
end
end

