function res = i_runsinglestep(t, N0, update)
%run a single step with the current odefile
global g_grind;
if update
   %update parameters!
   g_grind.hfun.curr = i_getodehandle(1,'');
end
if isempty(N0)
   res=[];
else
   siz1 = size(N0);
   switch length(siz1)
     case 2
       if g_grind.statevars.dim==1
           N0=reshape(N0,siz1(1) * siz1(2),1);
       end
     case 3
      N0 = reshape(N0, siz1(1) * siz1(2), siz1(3));
    case 4
      N0 = reshape(N0, siz1(1) * siz1(2) * siz1(3), siz1(4));
   end

   res = ones(size(N0));
   if numel(t)==1
       t=t+zeros(size(N0,1),1);
   elseif all(size(t)==siz1(1:2))
       t = reshape(t,size(N0,1),1);
   end

   if strcmp(g_grind.solver.opt.Vectorized,'on')
       yy = feval(g_grind.hfun.curr, transpose(t), transpose(N0));
       res = transpose(yy);
   else
       for i = 1:size(N0, 1)
          yy = feval(g_grind.hfun.curr, t(i), transpose(N0(i,:)));
          res(i, :) = transpose(yy);
       end

   end

   if length(siz1) > 2
      res = reshape(res, siz1);
   elseif length(siz1)== 2&&g_grind.statevars.dim==1&&numel(res)==prod(siz1)
      res=reshape(res,siz1);
   end

end

