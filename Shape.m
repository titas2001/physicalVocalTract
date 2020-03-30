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
%     [shapeType] 'linear' or 'exp' changes the type of curve after
%                  curveStartPos 
%                 'flat' creates a flat tube of amplitude minWidth. Ignores
%                 curveStartPos and maxWidth
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
        % sizes are in mm
    iArray = [12 12 32 32 32 32 32 32 24 16 10 10 10 12 14 24];
    eArray = [12 12 30 30 30 30 28 24 18 16 16 18 20 22 22 24];
    aArray = [12 12 26 16 12 14 20 26 30 34 38 38 34 30 28 32];
    oArray = [12 12 30 22 16 14 16 22 28 34 38 38 32 26 22 14];
    uArray = [12 12 30 30 30 30 26 18 14 22 26 24 22 20 14 16];
    
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
    %--------Vocal-tract-shapes-with-different-vocals---------------------------------
    if shapeType  == "i"
        leftover = N-floor(N/16)*16;
        for i = 1:16
            S(1+i*floor(N/16)-floor(N/16):i*floor(N/16)) = iArray(i)*0.001;
        end
        S(1+(i+1)*floor(N/16)-floor(N/16):i*floor(N/16)+leftover) = iArray(16)*0.001;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "e"
        leftover = N-floor(N/16)*16;
        for i = 1:16
            S(1+i*floor(N/16)-floor(N/16):i*floor(N/16)) = eArray(i)*0.001;
        end
        S(1+(i+1)*floor(N/16)-floor(N/16):i*floor(N/16)+leftover) = eArray(16)*0.001;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "a"
        leftover = N-floor(N/16)*16;
        for i = 1:16
            S(1+i*floor(N/16)-floor(N/16):i*floor(N/16)) = aArray(i)*0.001;
        end
        S(1+(i+1)*floor(N/16)-floor(N/16):i*floor(N/16)+leftover) = aArray(16)*0.001;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "o"
        leftover = N-floor(N/16)*16;
        for i = 1:16
            S(1+i*floor(N/16)-floor(N/16):i*floor(N/16)) = oArray(i)*0.001;
        end
        S(1+(i+1)*floor(N/16)-floor(N/16):i*floor(N/16)+leftover) = oArray(16)*0.001;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
    if shapeType  == "u"
        leftover = N-floor(N/16)*16;
        for i = 1:16
            S(1+i*floor(N/16)-floor(N/16):i*floor(N/16)) = uArray(i)*0.001;
        end
        S(1+(i+1)*floor(N/16)-floor(N/16):i*floor(N/16)+leftover) = uArray(16)*0.001;
        S = transpose(S);
        S = smoothdata(S, 'gaussian', 'SmoothingFactor', 0.05);
    end
end
