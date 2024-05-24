% Hann-square weights for linear fit at the edges, for MS1 smoothing.
function w = edgeWeights1(deg, m)
  beta = 0.65 + 0.35*exp(-0.55*(deg-4));
  fitLengthD = ((m+1)*beta)/(1+0.5*deg);
  fitLength = floor(fitLengthD);
  for i = 1 : fitLength+1
    cosine = cos(pi/2*(i-1)/fitLengthD);
    w(i) = cosine*cosine;
  end
end
