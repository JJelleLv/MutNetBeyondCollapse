%this function does an multivariate sensitivity analysis (Klepper 1998) on the  basis of
%the results of a sensitivity analysis. It calculates significance with dummy variables.
%There is no dependence on any GRIND function or variable
%
%mc_output= table with outputs (rows iterations columns model outputs (different variables/times))
%parsets = table with parameter values (rows iterations columns parameters)
%vlabels are labels for the model outputs (separately the output variable and time (nx2 cell matrix));
%plabels are parameter labels
%parndx ndx of parameters for dendrogram, empty=select by user
%
%
function [res,a,b]=i_mcsensanal(mc_output, parsets, vlabels, plabels,dendrpars,silent)
if nargout>0
    Soutput = sum_nonan(mc_output);
    Routput = max_nonan(mc_output) - min_nonan(mc_output);
    [res,a,b] = sensmatrix(parsets, mc_output, Soutput,Routput);
    return;
end

if nargin<5
    dendrpars={};
end

if nargin<6
    silent=false;
end

Ndummy = 1000;
Niter = size(mc_output, 1);
Nvars = size(mc_output, 2);
Npars = size(parsets, 2);
Soutput = sum_nonan(mc_output);
Routput = max_nonan(mc_output) - min_nonan(mc_output);
selectcases = Routput > 1E-12;
if ~all(selectcases)
   warning('GRIND:mcarlo:varsneglected','%d model outputs with negliglible variance removed\n', sum(~selectcases));
   f = find(selectcases);
   mc_output=mc_output(:,f);
   vlabels = vlabels(f, :);
   Nvars = size(mc_output, 2);
   Soutput = sum_nonan(mc_output);
   Routput = max_nonan(mc_output) - min_nonan(mc_output);
end

%sensitivity matrix
s_matrix = sensmatrix(parsets, mc_output, Soutput,Routput);

%dummy parameters are drawn at random between 0.9-1.1 (not sensitive to this choice).
dummy_matrix = sensmatrix(0.9 + rand(Niter, Ndummy) * 0.2, mc_output, Soutput,Routput);
plabels1 = cell(Npars, 1);
for j = 1:Npars
   plabels1{j} =  sprintf('%s', plabels{j, :});
end

vlabels1 = cell(Nvars, 1);
for j = 1:Nvars
   vlabels1{j} = sprintf('%s ', vlabels{j, :});
end

i_sensplot(mc_output,parsets,vlabels1,plabels1,silent);
i_klepperen(s_matrix,vlabels,plabels,dummy_matrix,'',dendrpars,silent);

% function [s_matrix,regress_a,regress_b] = sensmatrix(p, mc_output, Soutput,Routput)
% %create sensitivity matrix
% %Sp = sum(p);
% %Sp2 = sum(p.^2);
% m = sum(~isnan(mc_output),1);
% s_matrix = zeros(size(p,2), size(mc_output,2));
% if nargout>1
%     regress_a=s_matrix;
%     regress_b=s_matrix;
% end

% pRange = (max(p) - min(p))';
% for j = 1:size(p, 2)
%    %  pp = repmat(p(:, j), 1, size(mc_output,2)); %parameter values
%    % Sxy = sum(pp .* mc_output);
%    pp=p(:, j);
%    for i = 1: size(mc_output, 2)
%          Sxy = sum_nonan(p(:, j) .* mc_output(:, i));
%          Sp = sum(pp(~isnan(mc_output(:, i))));
%          Sp2 = sum(pp(~isnan(mc_output(:, i))).^2);
%          if m(i)>1
%             regress = [Sp2 Sp; Sp m(i)] \ [Sxy; Soutput(i)]; %linear regression (conform Grasman)
%          else
%              regress=[NaN,NaN];
%          end

%          a = regress(1);
%          if nargout>1
%              regress_a(j,i)=regress(1);
%              regress_b(j,i)=regress(2);
%          end

%          if Routput(i) == 0
%             s_matrix(j, i) = NaN; %delen door nul geeft nu NaN (anders zou het extreem gevoelig worden)
%          else
%             s_matrix(j, i) = a * pRange(j) / Routput(i);
%          end

%    end
% end


function [s_matrix,regress_a,regress_b] = sensmatrix(p, mc_output, Soutput,Routput)
%create sensitivity matrix
%Sp = sum(p);
%Sp2 = sum(p.^2);
m = sum(~isnan(mc_output),1);
s_matrix = zeros(size(p,2), size(mc_output,2));
regress_a=s_matrix;
regress_b=s_matrix;
pRange = transpose(max(p) - min(p));
for j = 1:size(p, 2)
   % vectorized for much more speed!
   pp=repmat(p(:,j),1,size(mc_output,2));
   Sxy = sum_nonan(repmat(p(:, j),1,size(mc_output, 2)) .* mc_output);
   Sp=sum(pp.*~isnan(mc_output));
   Sp2=sum((pp.*~isnan(mc_output)).^2);
  
   %Simple linear regression
   regress_a(j,:)=(Sxy-1./m.*Sp.*Soutput)./(Sp2-1./m.*Sp.^2);
   regress_b(j,:)=Soutput./m-regress_a(j,:).*Sp./m;
   s_mat = regress_a(j, :) .* pRange(j) ./ Routput;
   %replace zero output with NaN;
   s_mat(Routput == 0)=NaN;
   s_matrix(j, :)=s_mat;
end

function res = sum_nonan(A,p1)
if nargin==1
    p1=1;
end

A(isnan(A))=0;
res=sum(A,p1);

function res = max_nonan(A,p1,p2)
A(isnan(A))=-inf;
if nargin==1
  res=max(A);
elseif nargin==2
  res=max(A,p1);
elseif nargin==3
  res=max(A,p1,p2);
end


function res = min_nonan(A,p1,p2)
A(isnan(A))=inf;
if nargin==1
  res=min(A);
elseif nargin==2
  res=min(A,p1);
elseif nargin==3
  res=min(A,p1,p2);
end


