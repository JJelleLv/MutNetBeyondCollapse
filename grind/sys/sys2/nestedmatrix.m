%NESTEDMATRIX - create or evaluate nested network
%   Create a bipartite network as a matrix with a certain level of nestedness.
%   If you enter a low value for the nestedness, a random network is generated
%   intead. Algorithm by Jelle Lever and U.Bastolla (Bastolla et al 2009, Nature).
%
%   Usage:
%   NS = NESTEDMATRIX(A) - save the nestedness of matrix (A) to NS (algoritm
%      Bastolla).
%   A = NESTEDMATRIX(nrow,ncol,nestedness,connectivity) - create a nested
%      matrix with nrow rows and ncol columns, with a certain nestedness and
%      connectivity
%  [A, NS, conn] = NESTEDMATRIX(nrow,ncol,nestedness,connectivity,
%      forbiddenlinks) Create a nested network with forbidden links (see
%      Bastolla), save also the nestedness to NS and connectivity to conn)
function [A, NS,conn] = nestedmatrix(nrow,ncol,nestedness,connect,forbiddenlinks)
if (nargin==1)&&numel(nrow) > 1
    A1 = getnestedness(nrow);
    NS1 = sum(sum(nrow)) / numel(nrow);
    if nargout == 0
        fprintf('Size matrix:   [%dx%d]\n', size(nrow));
        fprintf('Nestedness:    %g\n', A1);
        fprintf('Connectedness: %g\n', NS1);
    else
        A = A1;
        NS = NS1;
    end

    return;
end

if nargin < 4
    connect = 0.5;
end

if nargin < 5
    forbiddenlinks = 0;
end

% Fraction of links before forbidden links are removed
CI_FL = connect / (1 - forbiddenlinks);
NS  = -1;
itrial = 1;
maxA = [];
maxNS = 0;
while NS < nestedness && itrial < 20
    %maximum 20 new networks if iter>1000
    NS = 999;
    i = 0;
    while NS > nestedness && i < 1000
        found = 0;
        while ~found
            % Make random matrix of mutualistic interactions
            A = rand(nrow, ncol);
            A=(A <= CI_FL);
            % Make matrix of forbidden links (zeros are forbidden)
            I_FL=rand(nrow, ncol) <= (1 - forbiddenlinks);
            
            %remove forbidden links
            A = A .* I_FL;
            
            %none of the species have no links
            found = all(sum(A, 1) > 0) & all(sum(A, 2) > 0);
        end
        NS = getnestedness(A);
        i = i + 1;
    end

    if NS > nestedness
        disp('nestedness too small for this algorithm, made a random matrix instead')
    end

    iter = 0;
    while (NS < nestedness) && (iter < (((nrow * ncol) * (nrow * ncol - 1)) / 2))
        A2 = A;
        % choose a species
        SP = ceil((nrow + ncol) * rand);
        ischanged = 0;
        if SP <= nrow % plant species
            i = SP;
            k = ceil(ncol * rand);
            l = ceil(ncol * rand);
            if A2(i, k)==1 && A2(i, l)==0 && sum(A2(:, k)) < sum(A2(:, l)) && I_FL(i, l)~=0
                A2(i, k) = A(i, l);
                A2(i, l) = A(i, k);
                ischanged = 1;
            elseif A2(i, k)==0 && A2(i, l)==1 && sum(A2(:, l)) < sum(A2(:, k)) && I_FL(i, k)~=0
                A2(i, k) = A(i, l);
                A2(i, l) = A(i, k);
                ischanged = 1;
            end
        else
            % animal species
            i = SP - nrow;
            k = ceil(nrow * rand);
            l = ceil(nrow * rand);
            if A2(k, i)==1 && A2(l, i)==0 && sum(A(k, :)) < sum(A(l, :)) && I_FL(l, i)~=0
                A2(k, i) = A(l, i);
                A2(l, i) = A(k, i);
                ischanged = 1;
            elseif A2(k, i)==0 && A2(l, i)==1 && sum(A2(l, :)) < sum(A2(k, :)) && I_FL(k, i)~=0
                A2(k, i) = A(l, i);
                A2(l, i) = A(k, i);
                ischanged = 1;
            end
        end
        if ischanged %this check makes the routine  > 10x faster
            % test that all species are still connected
            if ~(all(sum(A2, 1) > 0) && all(sum(A2, 2) > 0))
                %else discard
                A2 = A;
            end
            % calculate nestedness of new network
            NS = getnestedness(A2);
            A = A2;
        end

        iter = iter + 1;
    end
    itrial = itrial + 1;
    if NS > maxNS
        maxNS = NS;
        maxA = A;
    end

end

A = maxA;
NS = maxNS;
conn = sum(sum(A)) / (nrow * ncol);
% sort columns based on nr of interactions
[~, ndx] = sort(sum(A, 1), 'descend');
A = A(:, ndx);
[~, ndx] = sort(sum(A, 2), 'descend');
A = A(ndx, :);
disp(getnestedness(A));
%-------------------------------------------------------------------------%
% IMPLEMENTO EN MATLAB LA MEDIDA DE NESTEDNESS DEL ARTICULO DE BASTOLLA09 %
%-------------------------------------------------------------------------%
%
% code made by Ugo bastolla
% adapted by Egbert van Nes
%
%Se parte de una red de interacciones bipartita, animales en filas, y
%plantas en columnas.

function [NStotal,NSrows, NScols] = getnestedness(A)
%!!!!BE AWARE!!!! Plants and animals are switched
% normally: animales-horizontal, plantas vertical
%BINARIA=[1 1 1 1; 0 0 1 1; 0 0 1 1; 0 0 1 0; 0 0 1 0; 0 1 0 1; 0 0 0 1];
na = size(A, 1); %numero de animales
np = size(A, 2); %numero de plantas

sumNa = 0;
suma = 0;
%-SE CALCULA EL NESTEDNESS DE ANIMALES-%
for i = 1:na - 1 %se coge cada una de las filas y se ve el numero de unos que coinciden
    for j = i + 1:na
        coinciden = sum(A(j, :) .* A(i, :)); %nr of interactiones shared by two animales Bastolla: nij
        unafila = sum(A(j, :)); %es el numero de unos que tiene una de las filas escogidas Bastolla: ni
        otrafila = sum(A(i, :)); %es el numero de unos que tiene la otra de las filas escogidas Bastolla: nj
        
        % devide by min(ni,nj)
        if unafila<=otrafila && unafila~=0 && otrafila ~=0
            N = coinciden / unafila;
        elseif otrafila<=unafila && unafila~=0 && otrafila ~=0
            N = coinciden / otrafila;
        else
            N = 0;
        end%se divide el numero de unos en comun, por el valor minimo del numero de unos de las dos filas que estamos comparando
        sumNa = sumNa + N;
        suma = suma + 1;
    end

end
NSrows = sumNa / suma;
sumNp = 0;
sump = 0;
%save NESTEDNESSanimales
%-SE CALCULA EL NESTEDNESS DE PLANTAS-%
for i = 1:np - 1 %se coge cada una de las filas y se ve el numero de unos que coinciden
    for j = i + 1:np
        coinciden = sum(A(:, j) .* A(:, i)); %nr of interactiones shared by two animales Bastolla: nij
        unafila = sum(A(:, j)); %es el numero de unos que tiene una de las filas escogidas Bastolla: ni
        otrafila = sum(A(:, i)); %es el numero de unos que tiene la otra de las filas escogidas Bastolla: nj
        
        % devide by min(ni,nj)
        if unafila<=otrafila && unafila~=0 && otrafila ~=0
            N = coinciden / unafila;
        elseif otrafila<=unafila && unafila~=0 && otrafila ~=0
            N = coinciden / otrafila;
        else
            N = 0;
        end%se divide el numero de unos en comun, por el valor minimo del numero de unos de las dos filas que estamos comparando
        sumNp = sumNp + N;
        sump = sump + 1;
    end

end
NScols = sumNp / sump;
NStotal = (sumNa + sumNp) / (suma + sump);




