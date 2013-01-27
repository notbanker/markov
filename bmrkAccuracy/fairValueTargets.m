function trades = fairValueTargets(trades)
% function trades = fairValueTargets(trades)

nLeads = min(25,trades.nLeads(1));
nonPairedInterdealer = @(id_,isP_,rps_) deal (id_,false,4);
p = table.leading.valuesOfAnother(trades.trade_price,1:nLeads,nonPairedInterdealer,...
                                  trades.benchmark_id,trades.isPaired,trades.reporting_party_side);
s = table.leading.valuesOfAnother(trades.trade_size,1:nLeads,nonPairedInterdealer,...
                                  trades.benchmark_id,trades.isPaired,trades.reporting_party_side);                        
t = table.leading.valuesOfAnother(trades.received_time,1:nLeads,nonPairedInterdealer,...
                                  trades.benchmark_id,trades.isPaired,trades.reporting_party_side);                        
tRelative = t-repmat(trades.received_time,1,nLeads);

%% Simple fair value price
bj = cumsum(s,2)/1e6;
aj = [zeros(size(bj,1),1),bj(:,1:nLeads-1)];
w = s.*exp(-aj).*exp(-tRelative);
fvpNumer = nansum(p.*w,2);
fvpDenom = nansum(w,2);
trades.fair_value_price = fvpNumer./fvpDenom;
count = sum(~isnan(p),2);
trades.fair_value_price(count<3) = nan;

%% Slightly more consistent way
% We compute \int_m p(m) exp(-m) exp(-t(m)) dm where t(m), p(m) are
% piecewise constant in cumulative money m. We say [aj,bj] are the
% intervals in cumulative money. 
u = exp(-tRelative).*(exp(-aj)-exp(-bj));
consistentFvpNumer = nansum(p.*u,2);
consistentFvpDenom = nansum(u,2);
trades.fairer_value_price = consistentFvpNumer./consistentFvpDenom;
trades.fairer_value_price(count<3) = nan;


trades.fair_value_yield = approxYield(trades,'fair_value_price');
trades.fairer_value_yield = approxYield(trades,'fairer_value_price');



end