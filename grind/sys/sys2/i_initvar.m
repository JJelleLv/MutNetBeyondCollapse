% i_initvar = internal function to set initial values
function [N0, NP] = i_initvar(nowarn)
global g_grind;
 if g_grind.statevars.vector
   N0 = zeros(g_grind.statevars.dim,1);
   for j = 1:size(g_grind.statevars.vectnames, 2)
      avar = g_grind.statevars.vectnames{j};
      dim1 = g_grind.statevars.dims{j}.dim1;
      dim2 = g_grind.statevars.dims{j}.dim2;
      dim = dim1 * dim2;
      N01 = evalin('base', avar);
      if numel(N01) ~= dim
         if (nargin==0) || ~nowarn
            warning('GRIND:model:initialsize','The size of %s has been changed to match its size in the model definition\n', avar);
         end

         if numel(N01) == 1
            N01 = zeros(dim1, dim2) + N01;
         else
            N0a = zeros(dim1, dim2) + 0.001;
            if numel(N01) < dim
               N0a(1:size(N01, 1), 1:size(N01, 2)) = N01;
            else
               N0a = N01(1:dim1, 1:dim2);
            end

            N01 = N0a;
         end

         assignin('base', avar, N01);
      end

      N0(g_grind.statevars.dims{j}.from:g_grind.statevars.dims{j}.to) = N01(:);
   end

elseif ~isempty(g_grind.statevars.names)
    s=sprintf('%s;',g_grind.statevars.names{:});
    N0=evalin('base',sprintf('[%s];',s(1:end-1)));
    if length(N0)<g_grind.statevars.dim
        for i=1:g_grind.statevars.dim
            N1=evalin('base',g_grind.statevars.names{i});
            if isempty(N1)
                multassignin('base', g_grind.statevars.names{i},0.001);
            end

        end
        N0=evalin('base',sprintf('[%s];',s(1:end-1)));
    end
%    N0 = zeros(g_grind.statevars.dim, 1);
%    for i = 1:g_grind.statevars.dim
%       f=strfind(g_grind.statevars.names{i},'(');
%       if isempty(f)
%          N1 = evalin('base', g_grind.statevars.names{i});
%       else 
%          s = g_grind.statevars.names{i};
%          N1 = evalin('base',s(1:f(1)-1));
%          ndx=str2double(s(f(1)+1:strfind(g_grind.statevars.names{i},')')-1));
%          if length(N1)<ndx
%              N1=[];
%          end

%          if ~isempty(N1)
%             N1 = evalin('base', g_grind.statevars.names{i});
%          end; 
%       end

%       if ~isempty(N1)
%           N0(i) = N1(1);
%       else 
%          N0(i) = 0.001;
%          multassignin('base', g_grind.names{i},0.001);
%       end
%    end

else
    N0=[];
end

if nargout > 1
   if ~isempty(g_grind.permanent)
      NP = defpermanent('-p');
   else
      NP = [];
   end

end

