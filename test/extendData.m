% Extends the data by weighted linear extrapolation, for smoothing to the ends
function extData = extendData(data, m, fitWeights)
  datLength = length(data);
  fitLength = length(fitWeights);
  fitX = 1:fitLength;
  fitY = data(1:fitLength);
  [offset, slope] = fitWeighted(fitX, fitY, fitWeights);
  #fitBasis = [ones(1, fitLength); 1:fitLength];
  #[params] = LinearRegression (fitBasis', fitY', fitWeights');
  extData(1:m) = offset + (-m+1:0) * slope;
  extData(m+1 : datLength+m) = data;
  fitY = flip(data(datLength-fitLength+1 : datLength));
  [offset, slope] = fitWeighted(fitX, fitY, fitWeights);
  #[params] = LinearRegression (fitBasis', fitY', fitWeights');
  #extData(datLength+m+1 : datLength+2*m) = [ones(1, m); 0:-1:-m+1]' * params;
  extData(datLength+m+1 : datLength+2*m) = offset + (0:-1:-m+1) * slope;
end
