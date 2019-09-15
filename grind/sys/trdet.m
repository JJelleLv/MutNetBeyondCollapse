%TRDET   Plot of trace and determinant of the Jacobian and Hurwitz analysis
%  For 2D models plots of the trace and determinant are made. For higher dimensional models, an analysis of
%  the Routh-Hurwitz criteria is done. The results of this analysis are saved to a structure. 
%  The 2D trace and determinant of the Jacobian matrix determine the stability of equilibria. 
%  Plotting them help to determine the eigenvalues of <a href="trdet1.gif">difference</a> and 
%  <a href="trdet2.gif">differential equations</a>. This function works only for 2 dimensional equations.
%
%  Usage:
%  TRDET - derive trace and determinant plot.
%  TRDET TRUE - numerical Jacobian only
%  TRDET('jacobian',[1 2 3 4 1; 1 0 5 1 4; 1 2 3 4 5; 4 5 4 1 2 ;1 3 4 5 2]) - do the Hurwitz calculations for the provided Jacobian matrix.
%  res=TRDET - save the results of the analysis to a structure with fields:
%  res.F - the net feedback at level (x), see Levins, 1974 (should all be negative).
%  res.character - arguments of the characteristic equation
%  res.powers - powers of the characteristic equation
%  res.routh_table - the Routh table (number of changes in sign of the first column are the number of positive eigenvalues) 
%  res.hurwitz - Hurwitz determinants (should all be positive, to be stable) 
%  
%  TRDET('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'donumerical' [logical] - numerical or analytical Jacobian
%     'characteristic' [number] - enter values of the characteristic equation (for debugging)
%     'jacobian' [number] - values of the Jacobian
%  TRDET('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-d' - debug mode, some alternative algorithms are simulaneously calculated
%     
%         
%  See also eigen, findeq
%
%   Reference page in Help browser:
%      <a href="matlab:commands('trdet')">commands trdet</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function res1=trdet(varargin)
%(donumerical)
global g_grind;
fieldnams={'jacobian', 'n', 'values of the Jacobian',[];...
    'characteristic', 'n', 'enter values of the characteristic equation (for debugging)',[];...
    'donumerical', 'l', 'numerical or analytical Jacobian',false}';
args=i_parseargs(fieldnams,'donumerical','-d',varargin);
isdiffer=false;
N0=[];
if strcmp(args.opts,'-d') %debug mode
    debug=true;
else
    debug=false;
end
if isfield(args,'characteristic')&&~isempty(args.characteristic)
    args.jacobian=zeros(length(args.characteristic)-1);
else
    args.characteristic=[];
end
if ~isfield(args,'jacobian')||isempty(args.jacobian)
    i_parcheck;
    if ~isfield(args,'donumerical')
        args.donumerical = 0;
    end
    niters = g_grind.solver.iters;
    N0 = i_initvar;
    args.jacobian = i_eigen(args.donumerical,niters,N0);
    if g_grind.statevars.dim>2 && g_grind.solver.isdiffer
        error('grind:trdet:dim','trdet cannot deal with difference equations with more than 2 dimensions');
    end
    isdiffer=g_grind.solver.isdiffer;
end
J=args.jacobian;
if size(J,1)==2
    
    %     ud = get(get(0,'currentfigure'), 'userdata');
    %     if isfield(ud, 'iters')
    %         niters = ud.iters;
    %     end
    
    if any(any(isnan(J)))
        %    warning('grind:trdet:nanvals','Warning: Cannot determine trace and determinant as there are NaN values in the Jacobian');
        return;
    end
    if ~isempty(args.characteristic)
        tr=-args.characteristic(2);
        determ=args.characteristic(3);
    else
        tr = trace(J);
        determ = det(J);
        disp('Jacobian (J) =');
        disp(J);
    end
    fprintf('trace(J) = %g\n\n', tr);
    fprintf('det(J)   = %g\n', determ);
    hfig = i_makefig('trdet');
    hold off;
    plot(tr, determ, 'r+', 'markersize',10);
    i_plotdefaults(hfig);
    
    % daspect([1 1 1]);
    %  oldhold = ishold;
    hold on;
    set(hfig,'Name','Trace and determinant of Jacobian (J)');
    htitle=title('Trace and determinant of Jacobian (J)');
    set(htitle,'fontweight','normal');
    xlabel('trace(J)');
    ylabel('det(J)');
    H1 = get(hfig, 'CurrentAxes');
    max1 = max(abs([tr * 1.1 determ * 1.1]));
    if max1 < 0.18
        lim = [ -0.2 0.2];
    elseif max1 < 0.35
        lim = [ -0.4 0.4];
    elseif max1 < 0.45
        lim = [ -0.5 0.5];
    elseif max1 < 0.95
        lim = [ -1 1];
    elseif max1 < 1.4
        lim = [ -1.5, 1.5];
    elseif max1 < 10
        lim = [ -(round(max1) + 1), round(max1) + 1];
    else
        lim = [ -max1, max1];
    end
    set(H1, 'Xlim', lim);
    set(H1, 'Ylim', lim);
    x = lim(1):(lim(2) - lim(1)) / 100:lim(2);
    if isdiffer
        plot(x, 0.25 * x.^2,':k');
        x = [lim(1) lim(2)];
        plot(x, x - 1,':r');
        plot(x,  -x - 1,':r');
        x = [-2 2 0  -2];
        y = [1 1  -1 1];
        plot(x, y,'-k');
    else
        plot([lim(1),lim(2)],[0,0], 'k');
        plot([0,0],[lim(1),lim(2)], 'k');
        plot(x, 0.25 * x.^2,':k');
    end
    if ~isempty(N0)
        Nres1 = i_runsinglestep(1, N0', true)';
        if isdiffer
            Nres1 = Nres1 - N0;
        end
        if sum(Nres1) > 0.001
            warning('GRIND:trdet:noequilibrium','The model is not in an equilibrium\nUse <a href="matlab:findeq">findeq</a> to find the equilibrium');
        end
    end
end
if ~isdiffer
    %Routh-Hurwitz analysis
    
    %   Characteristic polynomial easy to check: filling the eigenvalues in the equation should be
    %   almost 0
    %
    dim=size(J,1);
    powers=(dim:-1:0);
    %cs=charpoly(J); %slower only available for MATLAB>2012 size(J) >40x40
    %cs=slowcharpoly(J); %really slow for size(J) > 20x20
    if isempty(args.characteristic)
        if debug
            try
                res.charpoly=charpoly(J);
            catch
                %if charpoly does not exist, just do nothing
            end
            if dim<20
               res.slowcharpoly=slowcharpoly(J); %very slow but exact
            end
            res.poly=poly(J); %very fast (uses eig) but less exact?
        end
        if dim>1500
            cs=poly(J); %very fast (uses eig) but less exact?
        else
            cs=fastcharpoly(J); %slow for size(J)>2000x2000
        end
        %cs=cs*cs(powers==dim);
    else
        cs=args.characteristic;
    end
    
    res.F=-cs(2:end);
    disp('Net feedback at each level (Levins)')
    for i=1:length(res.F)
        fprintf('F(%d) = %g\n',i,res.F(i));
    end
    disp('Characteristic equation:')
    s=sprintf('x^%d',dim);
    for i=2:length(powers)
        if cs(i)<0
            aplus='-';
            h=-1;
        else
            aplus='+';
            h=+1;
        end
        if powers(i)==1
            s=sprintf('%s %s %g*x',s,aplus,h*cs(i));
        elseif powers(i)==0
            s=sprintf('%s %s %g',s,aplus,h*cs(i));
        else
            s=sprintf('%s %s %g*x^%d',s,aplus,h*cs(i),powers(i));
        end
    end
    disp(s)
    [a,sp_c]=routh(cs);
    res.character=cs;
    res.powers=powers;
    res.routh_table=a;
    res.routh_special_cases=sp_c;
    %disp('Routh table:');
    %disp(a)
    chan=sum(abs(diff(sign(a(:,1)))/2));
    if chan==0
        fprintf('Stable: %d positive eigenvalues\n',chan);
    else
        fprintf('Unstable: %d positive eigenvalues\n',chan);
    end
    if debug
        disp('Hurwitz determinants:')
        
        res.Dx=hurwitz(cs);
        s=sprintf('Dx( 1)  =  %g\n',res.Dx(1));
        for i=2:length(res.Dx)
            s=sprintf('%sDx(%2d)  =  %g\n',s,i,res.Dx(i));
        end
        fprintf(s);
    end
    %Relationship between D and Routh table
    %A=D2/D1
    D=zeros(size(powers));
    D(1)=a(1);
    for i=2:length(D)
        D(i)=D(i-1)*a(i);
    end
    s=sprintf('D( 1)  =  %g\n',D(1));
    for i=2:length(D)
        s=sprintf('%sD(%2d)  =  %g\n',s,i,D(i));
    end
    fprintf(s);
    res.hurwitz=D;
    %error('GRIND:trdet:No2dim','This function can only be applied for two dimensional systems');
end
if nargout>0
    res1=res;
end
function res=hurwitz(cs)
%hurwitz matrices (without routh)
n=length(cs);
H=zeros(n);
cs1=cs;
if rem(n,2)==1
    cs1=[cs1 0];
    n_half=round((n+1)/2);
else
    n_half=round((n)/2);
end
cs_arr=flipud(reshape(cs1,[2,n_half]));
for i=1:2:n
    H(i:i+1,round(i/2):round(i/2)+n_half-1)=cs_arr;
end
res=zeros(1,n);
for i=1:length(res)
    res(i)=det(H(1:i-1,1:i-1));
end

function [rhTable,special_cases]=routh(coeffVector)
%% Routh-Hurwitz stability criterion
%
%  The Routh-Hurwitz stability criterion is a necessary (and frequently
%  sufficient) method to establish the stability of a single-input,
%  single-output(SISO), linear time invariant (LTI) control system.
%  More generally, given a polynomial, some calculations using only the
%  coefficients of that polynomial can lead us to the conclusion that it
%  is not stable.
%  Instructions
%  ------------
%
%  in this program you must give your system coefficients and the
%  Routh-Hurwitz table would be shown
%
%   Farzad Sagharchi ,Iran
%   2007/11/12
%   E van Nes, lot of bugs in special cases solved
%% Initialization
special_cases=struct('singlezero',false,'allzeros',false);
ceoffLength = length(coeffVector);
rhTableColumn = round(ceoffLength/2);
%  Initialize Routh-Hurwitz table with empty zero array
rhTable = zeros(ceoffLength,rhTableColumn);
%  Compute first row of the table
rhTable(1,:) = coeffVector(1,1:2:ceoffLength);
%  Check if length of coefficients vector is even or odd
if (rem(ceoffLength,2) ~= 0)
    % if odd, second row of table will be
    rhTable(2,1:rhTableColumn - 1) = coeffVector(1,2:2:ceoffLength);
else
    % if even, second row of table will be
    rhTable(2,:) = coeffVector(1,2:2:ceoffLength);
end
%% Calculate Routh-Hurwitz table's rows
%  Set epss as a small value
epss = 1E-8;
%  Calculate other elements of the table
for i = 3:ceoffLength
    %  special case 1: zero in the first column
    if (rhTable(i-1,1) == 0)&& ~all(rhTable(i-1,:) == 0)
        special_cases.singlezero=true;
        rhTable(i-1,1) = epss;
    end
    %  special case 2: row of all zeros
    if all(rhTable(i-1,:) == 0)
        special_cases.allzeros=true;
        order = ceoffLength - (i-2);
        cnt1 = 0;
        cnt2 = 1;
        RowBefore=rhTable(i-2,:);
        RowBefore=RowBefore./abs(RowBefore(1,1)); %normalize the first column
        for j = 1:rhTableColumn - 1
            rhTable(i-1,j) = (order - cnt1) * RowBefore(cnt2);
            cnt2 = cnt2 + 1;
            cnt1 = cnt1 + 2;
        end
    end

    %vectorized per row (speed bottleneck of routh)
    rhTable(i,1:end-1)=((rhTable(i-1,1) .* rhTable(i-2,2:end)) - ....
            (rhTable(i-2,1) .* rhTable(i-1,2:end))) ./ rhTable(i-1,1);
%
%     firstElemUpperRow = rhTable(i-1,1);
%     secondElemUpperRow = rhTable(i-2,1);
%
%     for j = 1:rhTableColumn - 1     
%         %  compute each element of the table
%         rhTable(i,j) = ((firstElemUpperRow * rhTable(i-2,j+1)) - ....
%             (secondElemUpperRow * rhTable(i-1,j+1))) / firstElemUpperRow;
%     end
        
end


function cs=slowcharpoly(J)
%Method for getting the charcteristic equation:
%    Bernard P. (2006) The coefficients of the characteristic polynomial in terms of the
%     eigenvalues and the elements of an n × n matrix. Applied Mathematics Letters 19 511–515        %very slow algorithm
dim=size(J,1);
powers=(dim:-1:0);
cs=zeros(size(powers));
for i=1:dim
    for j=1:length(powers)
        combs=nchoosek(1:dim,powers(j));
        if isempty(combs)
            cs(j)=det(J);
        else
            cs(j)=0;
            for k=1:size(combs,1)
                cs(j)=cs(j)+det_m_cols(J,combs(k,:));
            end
        end
    end
end
cs=cs*cs(powers==dim);

function res=det_m_cols(J,removecols)
J(removecols,:)=0;
J(:,removecols)=0;
for i=1:length(removecols)
    J(removecols(i),removecols(i))=-1;
end
res=det(J);

function p = fastcharpoly(A)
% Fast algorithm for characteristic polynomial of square matrices.
% Algortihm is described in
% La Budde's Method For Computing Characteristic Polynomials
% by Rizwana Rehman and Ilse C.F. Ipsen
%
% Input: square matrix A
%
% Output: charactersitic polynomial p of A
%
% Author: Sebastian J. Schlecht
% Date: 11.04.2015

%fastCharPoly - Fast algorithm for characteristic polynomial of square matrices.
%Algortihm is described in
%La Budde's Method For Computing Characteristic Polynomials
%by Rizwana Rehman and Ilse C.F. Ipsen
%
% Syntax:  p = fastCharPoly( A )
%
% Inputs:
%    A - square matrix complex or real
%
% Outputs:
%    p - coefficients of characteristic polynomial
%
% Example:
%    p = fastCharPoly( randn(4) )
%    p = charpoly(magic(4)) - fastCharPoly(magic(4))
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Dr.-Ing. Sebastian Jiro Schlecht,
% International Audio Laboratories, University of Erlangen-Nuremberg
% email address: sebastian.schlecht@audiolabs-erlangen.de
% Website: sebastianjiroschlecht.com
% 10. December 2018; Last revision: 10. December 2018


N = size(A, 1);
H = hess(A);
beta = diag(H(2:end,1:end-1));
alpha = diag(H);

pH = zeros(N+1,N+1);
pH(0+1, 1) = 1;
pH(1+1, 1:2) = [-alpha(1), 1];

for it = 2:N
    partB = flipud(cumprod(flipud(beta(1 : it-1))));
    partH = (H(1:it-1,it));
    partP = pH(1:it-1,:);
    
    rec = sum(bsxfun(@times, partB.*partH, partP),1);
    convp = conv(pH(it-1+1, :), [-alpha(it), 1]);
    
    pH(it+1, :) = convp(1:end-1) - rec;
end

p = fliplr(pH(end, :));
