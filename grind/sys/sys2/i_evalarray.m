%script i_evalarray
%g_l_s1 is the local variable with the commands
function g_res1=i_evalarray(g_l_s1)

%g_l_s1=parseeq(g_l_s1);
%g_l_s1=vectorize(g_l_s1);
g_l_obj=parsed_equation(g_l_s1);
g_l_s1=g_l_obj.vectorize;
clear('g_l_obj');
%g_l_s1=sprintf('%s',g_l_s1{:});
if nargout==0
    evalin('caller',g_l_s1);
else
   g_res1=evalin('caller',g_l_s1);
end

function g_res1=observed(ivar2)  %#ok<DEFNU>
global g_data;
    if ~isempty(g_data)
      if strcmp(ivar2,'t')
          g_res1=g_data.t;
      else 
          indx= strcmp(g_data.varlist,ivar2);
          g_res1 = g_data.obs(:,indx);
      end

    else
      g_res1=[nan;nan];
    end

% function s=vectorize(s) 
% ss=strcmp(s,'*');
% s(ss)={'.*'};
% ss=strcmp(s,'^');
% s(ss)={'.^'};
% ss=strcmp(s,'/');
% s(ss)={'./'};
% ss=strcmp(s,'&&');
% s(ss)={'&'};
% ss=strcmp(s,'||');
% s(ss)={'|'};
