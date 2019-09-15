function [Res,N] = i_changevar(S, Var1, Var2)
F = strfind(S, Var1);
lvar1 = length(Var1);
LF = sprintf('\n');  %#ok<SPRINTFN>
ops=[' .,$#@?=~&|><!{}[];:*^/\%+-()' LF]; %No convert
ops2=[ops ''''];
Res = S;
N=0;
for i = length(F):-1:1
   if ((F(i)==1) || ismember(S(F(i) - 1), ops)) && ...
         ((F(i)+lvar1 > length(S)) || ismember(S(F(i)+lvar1), ops2))
      Res = [S(1:F(i) - 1) Var2 Res(F(i) + lvar1:length(Res))];
      N=N+1;
   end

end


