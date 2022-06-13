function Tps = tpscore(Ws,axi)
    tnz = find(axi);
    if isequal(class(Ws),'cell')
        M=length(Ws);
        Tps = zeros(M,1);
        for m=1:M
            nz = find(Ws{m});
            FN = length(setdiff(tnz,nz));
            FP = length(setdiff(nz,tnz));
            TP = length(intersect(tnz,nz));
            Tps(m) = TP/(TP+FN+FP);
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
        Tps = zeros(M,1);
        for m=1:M
            if ~isempty(s)
                inds{s}=m;
            end
            nz = find(Ws(inds{:}));
            FN = length(setdiff(tnz,nz));
            FP = length(setdiff(nz,tnz));
            TP = length(intersect(tnz,nz));
            Tps(m) = TP/(TP+FN+FP);
        end
    else
        Tps=NaN;
    end
end
