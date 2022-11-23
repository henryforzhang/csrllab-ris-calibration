function y = wls_chan_est_fast_version(risEstEmplyed,Heff)
[~, oNb] = size(risEstEmplyed);
A = risEstEmplyed.';
const = 1/(1+oNb);
oneVec = ones(oNb,1);
t = A'*oneVec;
C1 = (A'*A-const*(t*t'))\A';
C2 = C1 - const*(C1*oneVec)*oneVec';
y = Heff*C2.';
end