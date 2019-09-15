% DRAWBINOMIAL - draw from a binomial distribution
% gets a value from a binomial distribution with probability P and number
% of independent experiments N. For for high numbers (N*P>10) or (N>1000000)
% it uses a Normal approximation.
%                              n     k      n-k
%         binomial =  P(k) = (   ) p   (1-p)
%                              k
% Usage:
% R = DRAWBINOMIAL(P,N) (P and N can also be a matrix)
function Result = drawbinomial(aP, N)
Result = zeros(size(N));
for i = 1:size(N, 1)
   for j = 1:size(N, 2)
      if length(aP) > 1
         if sum(size(aP)~=size(N))
            error('GRIND:drawbinomial:ArgError','Drawbinomial: size of matrices P and N should be the same');
         end
         P = aP(i, j);
      else
         P = aP;
      end
      if P > 0.5
         Result(i, j) = N(i, j) - drawbinomial(1-P, N(i, j));
      else
         aM = N(i, j) * P;
         if (aM > 10) || (N(i, j) >1000000)
            %large n -> centrale limietstelling}
            if aM * (1 - P) > 0
               Result(i, j) = round(aM + randn(1) * sqrt(aM * (1.0 - P)));
            else
               Result(i, j) = 0;
            end
         elseif aM <= 0
            Result(i, j) = 0;
         else
            %recursive notation of binomial:x);
            %           x=0    then P[x] = (1-p)^n
            %           x>=0   then P[x] = P[x-1]* (n-x+1)/x*p/(1-p)}
            Px = (1 - P)^N(i, j);
            CumP = Px;
            Uniform = rand(1);
            X = 0;
            PQ = P / (1 - P);
            while (CumP < Uniform) && (X<N(i, j))
               %may happen due to rounding errors
               X = X + 1;
               Px = Px * (N(i, j) - X + 1) / X * PQ;
               CumP = CumP + Px;
            end
            Result(i, j) = X;
         end
      end
   end
end
