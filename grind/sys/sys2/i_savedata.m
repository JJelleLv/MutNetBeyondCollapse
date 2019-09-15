function [varlist] = i_savedata(varlist, amatrix, checkonly)
global g_grind g_data;
if nargin < 3
   checkonly = 0;
end

hast = 0;
%err = 0;
datat=[];
ndx=strcmp(varlist,'###skip###')|strcmp(varlist,'');
for i = 1:length(varlist)
   if ~ndx(i)
      k = i_getno(varlist{i});
      if k.istime
         if hast
            %           err = 1;
            error('GRIND:setdata:TwoTimeCols','two columns represent time (t)');
         end

         hast = 1;
         ndx(i) = 1;
         datat = amatrix(:, i);
         if min(diff(datat))<0
               error('GRIND:setdata:DecrTime','Error in setdata: time must be monotonously increasing');         end

      elseif ~isempty(k.no) && k.ispar
         s = sprintf('Column %d (%s) represents a parameter, do you want to make it an external variable?' ...
            , i, varlist{i});
         button=questdlg(s,'data','Yes','Skip','Rename','Yes');
         switch button
          case 'Yes'
            answer=inputdlg({sprintf('What is the default value of %s?',varlist{i})},'data',1,{'0'});
            defextern(varlist{i}, answer);
          case 'Skip'
            ndx(i) = 1;
          case 'Rename'
            answer=inputdlg(sprintf('What is the new name of column %d?',i),'data',1,varlist(i));
            if ~isempty(answer)
               varlist{i} = answer;
            end

         end

      elseif isempty(k.no)
         if ~isvalidequation(varlist{i})
            s = sprintf('The variable or equation of column %d (%s) seems invalid, do you want to rename it?' ...
               ,i, varlist{i});
            button=questdlg(s,'data','Yes','Skip','Skip all unknown','Yes');
            if strcmp(button, 'Yes')
               answer=inputdlg({sprintf('What is the new name of column %d?',i)},'data',1,varlist(i));
               if ~isempty(answer)
                  varlist{i} = answer{1};
               end

            elseif strcmp(button, 'Skip all unknown')
               ndx(i) = 1;
               for j = i + 1:size(varlist, 2)
                  p = i_getno(varlist{j});
                  if isempty(p.no)&&~strcmp('t', varlist{j})&&~isvalidequation(varlist{j})
                     ndx(j) = 1;
                  end

               end

            else
               ndx(i) = 1;
            end

         end

      end

   end

end

nodata = 1;
for i = 1:length(varlist)
    if ~ndx(i)
        no = i_getno(varlist{i});
        if isempty(no.ndx)
            no.ndx=1;
        end

        if no.isext
            if isempty(datat)
                g_grind.externvars{no.no}.data = amatrix(:, i);
            elseif g_grind.externvars{no.no}.dim1*g_grind.externvars{no.no}.dim2==1
                g_grind.externvars{no.no}.data = [datat, amatrix(:,i)];
            else
                if size(g_grind.externvars{no.no}.data,1)~=length(datat) %we assume matrix external variables are in one matrix
                    g_grind.externvars{no.no}.data=nan+zeros(length(datat),g_grind.externvars{no.no}.dim1*g_grind.externvars{no.no}.dim2+1);
                end

                g_grind.externvars{no.no}.data(:,1) = datat;
                g_grind.externvars{no.no}.data(:,no.ndx+1) = amatrix(:,i);
            end

            ndx(i) = 1;
            g_grind.checks.lastsettings = {};
        elseif no.isvar
            nodata = 0;
        end

   end

end

for i=1:length(g_grind.externvars)
    ndx1=any(~isnan(g_grind.externvars{i}.data(:,2:end)), 2);
    g_grind.externvars{i}.data = g_grind.externvars{i}.data(ndx1, :);
end

varlist(ndx) = [];
amatrix(:, ndx) = [];
ndx2=any(~isnan(amatrix),2);
if ~isempty(datat)
  datat=datat(ndx2);
end

amatrix=amatrix(ndx2,:);
if ~checkonly&&~nodata
%    if ~isfield(g_data, 'options')
%       opts = optimset('fminsearch');
%       g_data.options.Display = opts.Display;
%       g_data.options.TolX = 1E-6;
%       g_data.options.TolFun = 1E-6;
%       g_data.options.MaxFunEvals1 = '200*numberOfVariables';
%       g_data.options.MaxIter1 = '200*numberOfVariables';
%       g_data.options.PosOnly = 0;
%       g_data.options.ResetEachStep = 0;
%    end

   g_data.varlist = varlist;
   g_data.obs = amatrix;
   g_data.pred = nan+zeros(size(amatrix));
   g_data.t = datat;
   g_data.minobs = min(g_data.obs);
   g_data.maxobs = max(g_data.obs);
end


function isok = isvalidequation(s)
try
   isok = 1;
   if ~strcontains(s,';')
       s=[s ';'];
   end

   evalin('base', s);
catch
   isok = 0;
end


