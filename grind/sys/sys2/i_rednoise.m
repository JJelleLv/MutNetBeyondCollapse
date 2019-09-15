function Redn=i_rednoise(N,labda,T0,beta)
Redn=zeros(N,1);
Redn(1)=T0;
for j=2:N
    Redn(j) = (1 - 1 / labda) * (Redn(j-1) - T0) + T0 + beta * randn;
end