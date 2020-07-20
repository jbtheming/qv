
function [QV, QVc, QVnoisy, Igamma] = QuadVar (data, data_type)
  ## Calculate realized quadratic variation of a series of prices
  N = length (data (:, 1)) - 1;
  Zraw = cell2mat (data (2:size (data, 1), 5));
  if strcmp (data_type, 'price')
    Z = log (Zraw);
  else
    Z = Zraw;
  endif
  ## tuning parameters
  k = 30;
  w = .47;
  lambda = 1 / 12;
  a0 = 12;
  s0 = std (Z (2:size (Z, 1)) - Z (1:size (Z, 1) - 1));
  del = min (datenum (data (3, 1)) - datenum (data (2, 1)), datenum (data (4, 1)) - datenum (data (3, 1)));
  u = k * del;
  v0 = a0 * s0 * (u ^ w);
  ## initial truncation parameter
  Rets_trunc = zeros (size (Z, 1) - 1, 1);
  for i = 1:size (Z, 1) - 1
    Rets_trunc (i, 1) = min(max(Z (i + 1, 1) - Z (i, 1), -v0), v0);
  endfor
  s = std (Rets_trunc);
  v = a0 * s * (u ^ w);
  theta = k * sqrt (del);
  ## create Zbar, locally averaged increments and Zhat, locally averaged
  ## squared increments
  Zbar = zeros (N - k + 1, 1);
  Zhat = zeros (N - k + 1, 1);
  for i = 1:N - k + 1
    for j = 1:k - 1
      Zbar (i, 1) = Zbar (i, 1) + g (j / k) * (Z (i + j) - Z (i + j - 1));
    endfor
    for j = 2:k
      Zhat (i, 1) = Zhat (i, 1) + ((g (j / k) - g ((j - 1) / k)) ^ 2) * ((Z (i + j - 1) - Z (i + j - 2)) ^ 2);
    endfor
  endfor
  ## calculate total quadratic variation QV, the continous part QVc and
  ## integrated noise variance Igamma
  QV = 0;
  QVc = 0;
  QVnoisy = 0;
  Igamma = 0;
  for i = 1:N - 1
    QVnoisy = QVnoisy + (Z (i + 1, 1) - Z (i, 1)) ^ 2;
  endfor
  for i = 1:N - k + 1
    QV = QV + Zbar (i) ^ 2 - (1 / 2) * Zhat (i);
    QVc = QVc + (Zbar (i) ^ 2 - (1 / 2) * Zhat (i)) * Ind (Zbar (i), v);
    Igamma = Igamma + Zhat (i) * Ind (Zbar (i), v);
  endfor
  QV = (1 / (k * lambda)) * QV;
  QVc = (1 / (k * lambda)) * QVc;
  Igamma = ((theta ^ 2) / (2 * k * lambda)) * Igamma;
  dummy=0
  ##
endfunction

