function [QCoV, QCoVc, corr, corr_c, corr_d] = QuadCoVar( data1, data2, data_type )
%Calculate realized quadratic covariation of 2 series of prices

%find common data sample window
t1 = datenum(data1(2,1));
T1 = datenum(data1(size(data1,1),1));
t2 = datenum(data2(2,1));
T2 = datenum(data2(size(data2,1),1));

t = max(t1,t2);
T = min(T1,T2);

data1 = data1(datenum(data1(2:end,1))<=T,:);
data1 = data1(datenum(data1(2:end,1))>=t,:);

data2 = data2(datenum(data2(2:end,1))<=T,:);
data2 = data2(datenum(data2(2:end,1))>=t,:);

N1 = length(data1(:,1)) - 1;
N2 = length(data2(:,1)) - 1;
N = min(N1,N2);
Z1raw = cell2mat(data1(2:size(data1,1),5));
Z2raw = cell2mat(data2(2:size(data2,1),5));
if strcmp(data_type, 'price')
    Z1 = log(Z1raw);
    Z2 = log(Z2raw);
else
    Z1 = Z1raw;
    Z2 = Z2raw;
end

%tuning parameters
k = 30;
w = .47;
lambda = 1/12;
a0 = 12; 
s10 = std(Z1(2:size(Z1,1)) - Z1(1:size(Z1,1)-1));
s20 = std(Z2(2:size(Z2,1)) - Z2(1:size(Z2,1)-1));
del = min(datenum(data1(3,1)) - datenum(data1(2,1)),datenum(data1(4,1)) - datenum(data1(3,1)));
u = k*del;
v10 = a0*s10*(u^w); %pre-truncation parameter
v20 = a0*s20*(u^w); %pre-truncation parameter
Rets_trunc1 = zeros(size(Z1,1) - 1, 1);
Rets_trunc2 = zeros(size(Z2,1) - 1, 1);
for i = 1:size(Z1,1) - 1
    Rets_trunc1(i,1) = min(Z1(i+1,1) - Z1(i,1),v10);
end
for i = 1:size(Z1,1) - 1
    Rets_trunc2(i,1) = min(Z2(i+1,1) - Z2(i,1),v20);
end
s1 = std(Rets_trunc1);
s2 = std(Rets_trunc2);
v1 = a0*s1*(u^w);
v2 = a0*s2*(u^w);
theta = k*sqrt(del);

%create Zbar, locally averaged increments and Zhat, locally averaged
%squared increments
Zbar1 = zeros(N - k + 1, 1);
Zbar2 = zeros(N - k + 1, 1);
Zhat = zeros(N - k + 1, 1);

for i = 1:N-k+1
   for j = 1:k-1
      Zbar1(i,1) = Zbar1(i,1) + g(j/k)*(Z1(i+j) - Z1(i+j-1));
      Zbar2(i,1) = Zbar2(i,1) + g(j/k)*(Z2(i+j) - Z2(i+j-1));
   end
   for j = 2:k
      Zhat(i,1) = Zhat(i,1) + ((g(j/k) - g((j-1)/k))^2)*((Z1(i+j-1) - Z1(i+j-2))*(Z2(i+j-1) - Z2(i+j-2)));
   end
end

%calculate total quadratic variation QV and the continous part QVc
QCoV = 0;
QCoVc = 0;

for i = 1:N-k+1
   QCoV = QCoV + (Zbar1(i)*Zbar2(i)) - (1/2)*Zhat(i);
   QCoVc = QCoVc + ((Zbar1(i)*Zbar2(i)) - (1/2)*Zhat(i))*Ind(Zbar1(i),v1)*Ind(Zbar2(i),v2);
end

QCoV = (1/(k*lambda))*QCoV;
QCoVc = (1/(k*lambda))*QCoVc;
QCoVd = QCoV - QCoVc;

[QV1,QVc1] = QuadVar(data1, data_type);
[QV2,QVc2] = QuadVar(data2, data_type);

corr = QCoV/(sqrt(QV1)*sqrt(QV2));
corr_c = QCoVc/(sqrt(QV1)*sqrt(QV2));
corr_d = QCoVd/(sqrt(QV1)*sqrt(QV2));

end


end

