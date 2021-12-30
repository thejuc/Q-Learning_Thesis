function NLL = SUM_GNG_loglikeli_action(Para,states,actions,rewards, Policy,QInit)

nSession = length(states);
nLL = nan(nSession,1);
for n=1:nSession
    nLL(n)  = GNG_loglikeli_action(Para,states{n},actions{n},rewards{n}, Policy,QInit);
end
NLL = sum(nLL,'omitnan');

end