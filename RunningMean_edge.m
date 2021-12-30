function [MTp]= RunningMean_edge(Tp,n)
%% running mean, average over n points, [mtp, std]
if isrow(Tp)
    Tp=Tp';
end
N=size(Tp,1);
mTp=nan(N,1);
stdTp=nan(N,1);
n2=floor(n/2);
if n<=N
    for i=1:N
        tp =nan;
        if i<=n2
            tp=Tp(1:i,1);
            
        end
        if i+n2>N
          %  tp=Tp(i:end,1);
            
        end
        if i>n2 && N-i>=n2
            tp=Tp(i-n2:i+n2,1);
        end
        stdtp=nanstd(tp);
        meantp=nanmean(tp);
       % tp=tp(find(tp>(meantp-5*stdtp) & tp<(meantp+5*stdtp))); % remove outliers
        mTp(i)=nanmean(tp);
        stdTp(i)=nanstd(tp);
        
    end
    
end

MTp=[mTp,stdTp];
end
