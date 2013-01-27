function trades = pairedTradeDetection(trades)
  
% Computes isPaired and adjusted sizes for paired trades
nLags = min(25,trades.nLags(1));

opposite = @(isD,bmk_id,size_cat) deal(~isD,bmk_id,size_cat);
timeOfPrevTwin = table.lagged.valuesOfAnother(trades.received_time,1,opposite,trades.isDealer,trades.benchmark_id,trades.size_category);
timeOfNextTwin = table.leading.valuesOfAnother(trades.received_time,1,opposite,trades.isDealer,trades.benchmark_id,trades.size_category);
same = @(isD,bmk_id,size_cat) deal(isD,bmk_id,size_cat);
timeOfPrevSame = table.lagged.valuesOfAnother(trades.received_time,1,same,trades.isDealer,trades.benchmark_id,trades.size_category);
timeOfNextSame = table.leading.valuesOfAnother(trades.received_time,1,same,trades.isDealer,trades.benchmark_id,trades.size_category);
isPairedAfter = timeOfPrevTwin>timeOfPrevSame & (trades.received_time-timeOfPrevTwin)<1;
isPairedBefore = timeOfNextTwin<timeOfNextSame & (timeOfNextTwin-trades.received_time)<15/(60*60*24);
trades.isPaired = isPairedAfter | isPairedBefore;
timeSincePrevTwin = trades.received_time-timeOfPrevTwin;
timeSincePrevTwin(~isPairedAfter) = 1000;
timeUntilNextTwin = timeOfNextTwin-trades.received_time;
timeUntilNextTwin(~isPairedBefore) = 1000;
tspt = fliplr(table.lagged.values(timeSincePrevTwin,1:nLags,trades.benchmark_id,trades.reporting_party_side));
tunt = fliplr(table.lagged.values(timeUntilNextTwin,1:nLags,trades.benchmark_id,trades.reporting_party_side));
useDD = @(id_,rps_) deal(id_,4);
sAdj = fliplr(table.lagged.valuesOfAnother(trades.trade_size, 1:nLags, useDD,trades.benchmark_id,trades.reporting_party_side));
pairDecayAfter = @(t_) ones(size(t_))-exp(-(40*t_+0.0001));
pairDecayBefore = @(t_) ones(size(t_))-exp(-(400*t_+0.0001));
SIZ = fliplr(table.lagged.valuesOfAnother(trades.trade_size, 1:nLags, useDD,trades.benchmark_id,trades.reporting_party_side));
sAdj(isPairedAfter,:) = SIZ(isPairedAfter,:).*pairDecayAfter(tspt(isPairedAfter,:));
sAdj(isPairedBefore,:) = SIZ(isPairedBefore,:).*pairDecayBefore(tunt(isPairedBefore,:));
trades.pair_adjusted_size = sAdj;

end