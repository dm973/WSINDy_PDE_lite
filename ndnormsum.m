function B = ndnormsum(varargin)
    B = varargin{1}.^2;
    for j=2:length(varargin)
        B = B+varargin{j}.^2;
    end
    B = sqrt(B);
end