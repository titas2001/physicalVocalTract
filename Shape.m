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
    %scaleFromCM2 = 1;
    scaleFromCM2 = 0.01;

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
    
    aJapanese2Array = [34 20 12 14 16 20 26 30 34 38 34 30 26 32];
    eJapanese2Array = [30 30 30 30 28 26 22 20 16 16 16 20 18 24];
    iJapanese2Array = [30 30 30 30 30 30 22 14 10 8 8 8 12 22];
    oJapanese2Array = [26 20 16 14 16 22 28 34 38 34 34 26 22 16];
    uJapanese2Array = [26 26 26 26 24 16 12 16 22 22 20 18 14 14];
    %~~~~~~~~~~~~~~~~Japanese~Vowels~in~mm~~~~~~~~~~~~~~~~~~~~~~
    iJapaneseArray = [12 12 32 32 32 32 32 32 24 16 10 10 10 12 14 24];
    eJapaneseArray = [12 12 30 30 30 30 28 24 18 16 16 18 20 22 22 24];
    aJapaneseArray = [12 12 26 16 12 14 20 26 30 34 38 38 34 30 28 32];
    oJapaneseArray = [12 12 30 22 16 14 16 22 28 34 38 38 32 26 22 14];
    uJapaneseArray = [12 12 30 30 30 30 26 18 14 22 26 24 22 20 14 16];
    %calculate area funcion in cm2
    iJapanese2Array = (iJapaneseArray/20).^2*pi;
    eJapanese2Array = (eJapaneseArray/20).^2*pi;
    aJapanese2Array = (aJapaneseArray/20).^2*pi;
    oJapanese2Array = (oJapaneseArray/20).^2*pi;
    uJapanese2Array = (uJapaneseArray/20).^2*pi;
    %~~~~~~~~~~~~~~~~American~Vowels~in~cm2~~~~~~~~~~~~~~~~~~~~~~
    aAmericanArray  = [1.56 3.10 3.74 2.48 1.28 0.60 0.73 1.28 1.39 1.31 1.43 1.90 3.22 4.44 4.83 3.89 4.72 2.03 2.49 2.82];
    aeAmericanArray = [1.91 2.69 4.53 2.05 1.33 1.44 2.82 4.54 4.04 3.24 2.82 3.20 4.32 5.13 4.17 4.98 6.31 5.65 7.03];
    iAmericanArray  = [1.20 2.73 4.11 4.63 6.05 7.62 7.64 7.99 7.09 4.55 2.63 1.81 1.10 0.69 0.85 0.80 0.50 1.10 1.65];
    uAmericanArray  = [3.11 5.18 6.44 6.05 5.76 6.20 5.19 3.72 3.06 2.31 1.05 1.00 0.67 0.92 1.57 2.30 3.78 4.24 3.89 2.13 0.66];
    %~~~~~~~~~~~~~~~~British~Vowels~in~cm2~~~~~~~~~~~~~~~~~~~~~~
    aBritishArray  = [2.02 3.34 2.56 1.18 0.76 0.67 0.81 0.80 1.47 2.48 2.85 2.76 3.29 3.60 3.10 2.90 2.53 4.23];
    aeBritishArray = [0.48 1.86 1.98 1.92 0.83 0.85 1.68 1.56 1.64 2.34 2.66 2.50 1.77 2.19 2.28 2.35 4.82 7.89];
    iBritishArray  = [1.33 1.87 4.31 5.38 6.36 8.11 7.73 6.77 5.68 4.35 2.93 1.64 1.01 0.55 0.54 0.56 1.62 2.37 3.33];
    uBritishArray  = [0.55 1.92 4.52 6.87 6.87 7.12 6.08 5.09 4.48 2.80 2.11 1.53 0.74 1.20 0.79 1.30 2.03 2.79 3.36 2.44 1.07];
    % ref 
    % Analysis of vocal tract shape and dimensions using magnetic
    % resonance imaging: Vowels 
    % by T. Baer,J.C. Gore, L.C. Gracco and P.W. Nye 
    %~~~~~~~~~~~~~~~~French~Vowels~in~cm2~~~~~~~~~~~~~~~~~~~~~~
    aFrenchArray = [1.8 1.8 2.8 1.5 0.8 1.3 1.5 1.7 2.8 4.5 7.1 9.3 13.5 15.5 7.8];
    iFrenchArray = [3.0 3.7 4.5 6.6 9.3 10.3 9.4 7.5 4.5 1.4 0.6 0.4 0.3 2.2];
    uFrenchArray = [1.4 3.4 6.6 8.7 9.8 9.2 7.6 5.8 4.4 2.3 0.7 1.5 3.9 8.6 6.9 1.5 1.3];
    % ref
    % Vocal Tract Area Function for Vowels Using
    % Three-Dimensional Magnetic Resonance
    % Imaging. A Preliminary Study
    % by Philippe Cle´ment, Ste´phane Hans, Dana M. Hartl, Shinji Maeda, 
    % Jacqueline Vaissie`re, and Daniel Brasnu
    
    
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
    if shapeType  == "aJ" % a Japanese
        L = length(aJapanese2Array);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = aJapanese2Array(i) * scaleFromCM2;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = aJapanese2Array(L) * scaleFromCM2;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "eJ"
        L = length(eJapanese2Array);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = eJapanese2Array(i) * scaleFromCM2;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = eJapanese2Array(L) * scaleFromCM2;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "iJ"
        L = length(iJapanese2Array);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = iJapanese2Array(i) * scaleFromCM2;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = iJapanese2Array(L) * scaleFromCM2;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "oJ"
        L = length(oJapanese2Array);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = oJapanese2Array(i) * scaleFromCM2;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = oJapanese2Array(L) * scaleFromCM2;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "uJ"
        L = length(uJapanese2Array);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = uJapanese2Array(i) * scaleFromCM2;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = uJapanese2Array(L) * scaleFromCM2;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %--------Vocal-tract-shapes-with-different-American-vocals---------------------------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if shapeType  == "aA" % a American
        L = length(aAmericanArray);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = aAmericanArray(i) * scaleFromCM2;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = aAmericanArray(L) * scaleFromCM2;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "aeA"
        L = length(aeAmericanArray);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = aeAmericanArray(i) * scaleFromCM2;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = aeAmericanArray(L) * scaleFromCM2;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "iA"
        L = length(iAmericanArray);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = iAmericanArray(i) * scaleFromCM2;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = iAmericanArray(L) * scaleFromCM2;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "uA"
        L = length(uAmericanArray);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = uAmericanArray(i) * scaleFromCM2;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = uAmericanArray(L) * scaleFromCM2;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %--------Vocal-tract-shapes-with-different-British-vocals----------------------------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if shapeType  == "aB" % a British
        L = length(aBritishArray);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = aBritishArray(i) * scaleFromCM2;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = aBritishArray(L) * scaleFromCM2;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "aeB"
        L = length(aeBritishArray);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = aeBritishArray(i) * scaleFromCM2;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = aeBritishArray(L) * scaleFromCM2;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "iB"
        L = length(iBritishArray);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = iBritishArray(i) * scaleFromCM2;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = iBritishArray(L) * scaleFromCM2;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "uB"
        L = length(uBritishArray);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = uBritishArray(i) * scaleFromCM2;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = uBritishArray(L) * scaleFromCM2;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------Vocal-tract-shapes-with-different-British-vocals----------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if shapeType  == "aF" % a French
        L = length(aFrenchArray);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = aFrenchArray(i) * scaleFromCM2;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = aFrenchArray(L) * scaleFromCM2;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "iF"
        L = length(iFrenchArray);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = iFrenchArray(i) * scaleFromCM2;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = iFrenchArray(L) * scaleFromCM2;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "uF"
        L = length(uFrenchArray);
        leftover = N-floor(N/L)*L;
        for i = 1:L
            S(1+i*floor(N/L)-floor(N/L):i*floor(N/L)) = uFrenchArray(i) * scaleFromCM2;
        end
        S(1+(i+1)*floor(N/L)-floor(N/L):i*floor(N/L)+leftover) = uFrenchArray(L) * scaleFromCM2;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
end
