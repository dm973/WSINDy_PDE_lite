function dW = wnorm(Ws,axi,p)
    if ~isempty(axi)
        eq=size(axi,2);
        if isequal(class(Ws),'cell')
            M=length(Ws);
            dW = zeros(M,eq);
            for m=1:M
                W = Ws{m};
                dW(m,:) = vecnorm(W-axi,p,1)./vecnorm(axi,p,1);
            end
        elseif isequal(class(Ws),'double')
            dimsW = size(Ws);
            dimsA = size(axi);
            s = find(~ismember(dimsW,dimsA));
            inds = repmat({':'},1,length(dimsW));
            if isempty(s)
                M=1;
            else
                M=dimsW(s);
            end            
            dW = zeros(M,eq);
            for m=1:M
                if ~isempty(s)
                    inds{s}=m;
                end
                W = Ws(inds{:});
                dW(m,:) = vecnorm(W-axi,p,1)./vecnorm(axi,p,1);
            end
        end
    else
        dW=[];
    end
end
