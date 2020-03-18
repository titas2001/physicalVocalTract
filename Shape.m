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
%     Shape(N,0.5, 4, 'linear') will get function with linear curve with max
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
end

