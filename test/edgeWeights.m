% Hann-square weights for linear fit at the edges, for MS smoothing.
function w = edgeWeights(deg, m)
  beta = 0.70 + 0.14*exp(-0.6*(deg-4));
  fitLengthD = ((m+1)*beta)/(1.5+0.5*deg);
  fitLength = floor(fitLengthD);
  for i = 1 : fitLength+1
    cosine = cos(pi/2*(i-1)/fitLengthD);
    w(i) = cosine*cosine;
  end
end
