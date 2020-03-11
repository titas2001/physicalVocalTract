function [S] = Shape(N, curveStartPos,  maxWidth, shapeType )
 
%
%     Shape function inputs:
%     [N] determines the length of the tube
% 
%     [curveStartPos] value betweem 0 and 1, determines where selected curve starts 
%         in function 0 is start and 1 is end
% 
%     [maxWidth] determines width of the tube in the open end of the tube 
% 
%     [shapeType] 'linear' or 'exp' changes the type of curve after curveStartPos 
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
        S = ones(N,1);
        D = ceil(N*curveStartPos);
        expPower = (log(maxWidth))/(log(exp(1))*(N-D));
        for i = D:1:N
            S(i) = exp((i-D)*expPower);
        end
    end
    
    % linear curve shape
    if shapeType  == "linear"
        S = ones(N,1);
        D = floor(N*curveStartPos);
        c = (maxWidth-1)/(N-D);
        for i = D:N
            S(i) = (i-D)*c + 1;
        end
    end
end

