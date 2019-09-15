%implements the method of O Klepper to do a cluster analysis on a
%sensitivity matrix
function i_klepperen(s_matrix, vlabels, plabels, dummy_matrix, atitle, dendrpars, silent)
if nargin < 5
   atitle = '';
end

if nargin<6
    dendrpars={};
end

if nargin<7
    silent=false;
end

Npars = size(plabels, 1);
Nvars = size(vlabels, 1);
Ndummy =  size(dummy_matrix, 1);
nvlabel =  size(vlabels, 2); %you may have more than one variable label (e.g. var and time)
nplabel =  size(plabels, 2); %you may have more than one parameter label
if exist('i_figno','file')
   nfig = i_figno('mcarlo');
else
   nfig = 1;
end

if ~isempty(dummy_matrix)
   Lidummy = sqrt(sum(dummy_matrix.^2, 2));
   
   hfig=i_figure(nfig);
   if silent
       set(hfig,'visible','off');
   end
   if any(~isreal(Lidummy(:)))
       Lidummy=real(Lidummy);
   end
   hist(Lidummy, 30);
   i_plotdefaults(hfig);
   htitle=title('Distribution total sensitivity of dummy parameters');
   set(htitle,'fontweight','normal');
   xlabel('total sensitivity');
   ylabel('frequency');
   drawnow;
   percentiles = [0.01 0.05 0.1 0.5 0.9, 0.95, 0.99];
   signif = i_makepercentiles(transpose(Lidummy), percentiles);
   percdummy = transpose(i_makepercentiles(transpose(dummy_matrix), percentiles));
   fprintf('Analysis of %d dummy parameters\n', Ndummy);
   for i = 1:length(percentiles)
      if percentiles(i) > 0.5
         fprintf('%g%% percentile: %g\n', percentiles(i) * 100, signif(i));
      end

   end

else
   percentiles = [];
   percdummy = [];
   signif = [];
end

% Generate a nice table and copy to clipboard
if ~isempty(percentiles)
   SensitivityMatrix = cell(Npars + 1+nvlabel+length(percentiles), Nvars + 1+nplabel);
else
   SensitivityMatrix = cell(Npars + nvlabel, Nvars + 1+nplabel);
end
% (2) analyse the sensitivity matrix (statistical toolbox needed?)

Li = sqrt(sum_nonan(s_matrix.^2, 2));
plabels1 = cell(Npars, 1);
plabels2 = plabels1;

for j = 1:Npars
   for i = 1:nplabel
      SensitivityMatrix{j + nvlabel, 1+i} = plabels{j,i};
   end

   SensitivityMatrix{j +  nvlabel, 1} = Li(j);
   for i = 1:Nvars
       %vectorization does not help for performance as num2cell is not
       %vectorized
      SensitivityMatrix{j + nvlabel, i + nplabel+1} = s_matrix(j,i);
   end

   plabels1{j} =  sprintf('%s', plabels{j, :});
   plabels2{j} = sprintf('%5.4f    %s', Li(j), plabels1{j,:});
end

vlabels1 = cell(Nvars, 1);
for j = 1:nvlabel
   for i = 1:Nvars
      SensitivityMatrix{j, i + nplabel+1} = vlabels{i,j};
      vlabels1{i} = sprintf('%s ', vlabels{i, :});
   end

end

for i = 1:length(percentiles)
   SensitivityMatrix{Npars + 1+nvlabel+i, 2} = sprintf('%g%% percentile', percentiles(i) * 100);
   SensitivityMatrix{Npars + 1+nvlabel+i, 1} = signif(i);
   for j = 1:Nvars
      SensitivityMatrix{Npars+ 1+nvlabel+i, j + nplabel+1} = percdummy(i,j);
   end

end


if exist('uitable','builtin')
    if ~silent
        if isempty(atitle)
            atitle = 'Full sensitivity matrix';
        end

        i_table(SensitivityMatrix, atitle);
    end
else
    varcopy(SensitivityMatrix); %this is the GRIND clipboard function
    disp('The Sensitivity Matrix is copied to the clipboard, paste for instance in Excel')
end

if ~silent
  i_sensmatrixmovie(s_matrix, vlabels, plabels);
  i_senstimeplot( s_matrix,vlabels,plabels1,percdummy);
end

%plots of each output versus each parameter
%make dendrogram
if ~i_hastoolbox('stats')&&length(dendrpars)~=1
   error('GRIND:mcarlo:NoStatToolbox','mcarlo: statistics toolbox required for the cluster analysis');
end

if length(plabels2) > 2
   if isempty(dendrpars)
     parndx = i_parseldlg( plabels2, Li, percentiles, signif);
   else
      parndx=find(ismember(plabels1,dendrpars));
   end

   %nonans=all(~isnan(scaled));
   %scaled = scaled(parndx, nonans);
   nonans = all(~isnan(s_matrix));
   len = 0;
   plabels1 = plabels1(parndx);
   for i = 1:length(plabels1)
      if len < length(plabels1{i})
         len = length(plabels1{i});
      end

   end

   plabels2 = plabels2(parndx);
   Li = Li(parndx);
   for i = 1:length(plabels1)
      if len < length(plabels1{i})
         len = length(plabels1{i});
      end

      plabels2{i}=sprintf('%5.4f %s%s', Li(i),repmat(' ',1,len+1-length(plabels1{i})),plabels1{i});
   end

   s_matrix = s_matrix(parndx, nonans);
   nlab = length(plabels2);
   if nlab < 2
      warning('GRIND:mcarlo:dendrogram1var','Cannot make a dendrogram if there is only one parameter/variable');
   else
      Y = pdist(s_matrix, str2func('sinedist')); % calculate distances
      if ~all(isnan(Y))
          Z = linkage(Y, 'average');
          hfig=i_figure(nfig + 3);
          if silent
              set(hfig,'visible','off');
          end

          dendrogram(Z,nlab,'orientation','right','labels',plabels2);
          %i_plotdefaults(hfig);
          xlabel('Sine distance');
          set(gca, 'fontname','Courier New')
      end

   end

end


function [d] = sinedist(u, V)
%absolute sine is sqrt(1-cos)
sumX2 = sum_nonan(u.^2, 2);
sumY2 = sum_nonan(V.^2, 2);
sumXY = sum_nonan(V .* repmat(u, size(V, 1), 1), 2);
d = sqrt(1 - (sumXY ./ sqrt(sumX2 .* sumY2)).^2);

function res = sum_nonan(A, p1)
if nargin == 1
   p1 = 1;
end
A(isnan(A)) = 0;
res = sum(A, p1);



