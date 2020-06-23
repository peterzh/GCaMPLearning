function out = grammCallbackBinomialConfidenceInterval(y)
% This function is used by gramm to generate binomial confidence intervals
% for psychometric curves. 
% Peter Zatka-Haas
[p, pci] = binofit( sum(y), length(y));
out = [p; pci'];
end