function [S] = Shape(N, curveStartPos, minWidth, maxWidth, shapeType)
 
%
%     Shape function inputs:
%     [N] determines the length of the tube
% 
%     [curveStartPos] value betweem 0 and 1, determines where selected curve starts 
%         in function 0 is start and 1 is end
%       
%     [minWidth] determines width of the tube in the flat part
%
%     [maxWidth] determines width of the tube in the open end of the tube 
% 
%     [shapeType] - 'linear' or 'exp' changes the type of curve after
%                    curveStartPos 
%                 - 'flat' creates a flat tube of amplitude minWidth. Ignores
%                   curveStartPos and maxWidth
%                 - 'a', 'e', 'i', 'o', 'u' vowels shapes. Ignores 
%                   curveStartPos
% 
%     example:
% 
%     Shape(N,0.5,0.5, 4, 'linear') will get function with linear curve with max
%     value 4 and the curve will start at halfway point of N
%
    % check if inputs are valid
    if floor(N) ~= ceil(N) || N<0
        error('N has to be a positive integer');
    end
    if curveStartPos <= 0 || curveStartPos >= 1
        error('curveStatPos has to be a double between 0 and 1');
    end
    if maxWidth <= 0
        error('maxWidth has to be a possitive number');
    end
    %~~~~~~~~~~~~~~~~Japanese~Vowels~in~mm~~~~~~~~~~~~~~~~~~~~~~
    iArray = [12 12 32 32 32 32 32 32 24 16 10 10 10 12 14 24];
    eArray = [12 12 30 30 30 30 28 24 18 16 16 18 20 22 22 24];
    aArray = [12 12 26 16 12 14 20 26 30 34 38 38 34 30 28 32];
    oArray = [12 12 30 22 16 14 16 22 28 34 38 38 32 26 22 14];
    uArray = [12 12 30 30 30 30 26 18 14 22 26 24 22 20 14 16];
    
    %~~~~~~~~~~~~~~~~American~Vowels~in~cm~~~~~~~~~~~~~~~~~~~~~~
    aAmericanArray  = [1.56 3.10 3.74 2.48 1.28 0.60 0.73 1.28 1.39 1.31 1.43 1.90 3.22 4.44 4.83 3.89 4.72 2.03 2.49 2.82];
    aeAmericanArray = [1.91 2.69 4.53 2.05 1.33 1.44 2.82 4.54 4.04 3.24 2.82 3.20 4.32 5.13 4.17 4.98 6.31 5.65 7.03];
    iAmericanArray  = [1.20 2.73 4.11 4.63 6.05 7.62 7.64 7.99 7.09 4.55 2.63 1.81 1.10 0.69 0.85 0.80 0.50 1.10 1.65];
    uAmericanArray  = [3.11 5.18 6.44 6.05 5.76 6.20 5.19 3.72 3.06 2.31 1.05 1.00 0.67 0.92 1.57 2.30 3.78 4.24 3.89 2.13 0.66];
    %~~~~~~~~~~~~~~~~British~Vowels~in~cm~~~~~~~~~~~~~~~~~~~~~~
    aBritishArray  = [2.02 3.34 2.56 1.18 0.76 0.67 0.81 0.80 1.47 2.48 2.85 2.76 3.29 3.60 3.10 2.90 2.53 4.23];
    aeBritishArray = [0.48 1.86 1.98 1.92 0.83 0.85 1.68 1.56 1.64 2.34 2.66 2.50 1.77 2.19 2.28 2.35 4.82 7.89];
    iBritishArray  = [1.33 1.87 4.31 5.38 6.36 8.11 7.73 6.77 5.68 4.35 2.93 1.64 1.01 0.55 0.54 0.56 1.62 2.37 3.33];
    uBritishArray  = [0.55 1.92 4.52 6.87 6.87 7.12 6.08 5.09 4.48 2.80 2.11 1.53 0.74 1.20 0.79 1.30 2.03 2.79 3.36 2.44 1.07];
    % ref 
    % Analysis of vocal tract shape and dimensions using magnetic
    % resonance imaging: Vowels 
    % by T. Baer,J.C. Gore, L.C. Gracco and P.W. Nye 
    % exponential curve shape
    if  shapeType == "exp" 
        S = ones(N,1) * minWidth;
        D = ceil(N*curveStartPos);
        expPower = log(1-minWidth+maxWidth)/(N-D-1);
        for i = 1:(N-D)
            S(i+D) = exp((i-1)*expPower) - 1 + minWidth;
        end
    end
    
    % linear curve shape
    if shapeType  == "linear"
        S = ones(N,1) * minWidth;
        D = floor(N*curveStartPos);
        c = (maxWidth-minWidth)/(N-D);
        for i = D:N
            S(i) = (i-D)*c + minWidth;
        end
    end
    
    if shapeType == "flat"
        S = ones(1,N) * minWidth;
    end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %--------Vocal-tract-shapes-with-different-Japanese-vocals---------------------------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if shapeType  == "iJ"
        leftover = N-floor(N/16)*16;
        for i = 1:16
            S(1+i*floor(N/16)-floor(N/16):i*floor(N/16)) = iArray(i)*0.001;
        end
        S(1+(i+1)*floor(N/16)-floor(N/16):i*floor(N/16)+leftover) = iArray(16)*0.001;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "eJ"
        leftover = N-floor(N/16)*16;
        for i = 1:16
            S(1+i*floor(N/16)-floor(N/16):i*floor(N/16)) = eArray(i)*0.001;
        end
        S(1+(i+1)*floor(N/16)-floor(N/16):i*floor(N/16)+leftover) = eArray(16)*0.001;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "aJ"
        leftover = N-floor(N/16)*16;
        for i = 1:16
            S(1+i*floor(N/16)-floor(N/16):i*floor(N/16)) = aArray(i)*0.001;
        end
        S(1+(i+1)*floor(N/16)-floor(N/16):i*floor(N/16)+leftover) = aArray(16)*0.001;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "oJ"
        leftover = N-floor(N/16)*16;
        for i = 1:16
            S(1+i*floor(N/16)-floor(N/16):i*floor(N/16)) = oArray(i)*0.001;
        end
        S(1+(i+1)*floor(N/16)-floor(N/16):i*floor(N/16)+leftover) = oArray(16)*0.001;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "uJ"
        leftover = N-floor(N/16)*16;
        for i = 1:16
            S(1+i*floor(N/16)-floor(N/16):i*floor(N/16)) = uArray(i)*0.001;
        end
        S(1+(i+1)*floor(N/16)-floor(N/16):i*floor(N/16)+leftover) = uArray(16)*0.001;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %--------Vocal-tract-shapes-with-different-American-vocals---------------------------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if shapeType  == "aA" % a american
        L = length(aAmericanArray);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = aAmericanArray(i)*0.01;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = aAmericanArray(L)*0.01;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "aeA"
        L = length(aeAmericanArray);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = aeAmericanArray(i)*0.01;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = aeAmericanArray(L)*0.01;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "iA"
        L = length(iAmericanArray);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = iAmericanArray(i)*0.01;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = iAmericanArray(L)*0.01;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "uA"
        L = length(uAmericanArray);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = uAmericanArray(i)*0.01;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = uAmericanArray(L)*0.01;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %--------Vocal-tract-shapes-with-different-British-vocals----------------------------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if shapeType  == "aB" % a american
        L = length(aBritishArray);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = aBritishArray(i)*0.01;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = aBritishArray(L)*0.01;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "aeB"
        L = length(aeBritishArray);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = aeBritishArray(i)*0.01;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = aeBritishArray(L)*0.01;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "iB"
        L = length(iBritishArray);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = iBritishArray(i)*0.01;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = iBritishArray(L)*0.01;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "uB"
        L = length(uBritishArray);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = uBritishArray(i)*0.01;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = uBritishArray(L)*0.01;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
end
