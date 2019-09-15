function [sizes1, freq1, B1] = powerlaw(A, critA)
if nargin < 2
   critA = mean(mean(A));
end
C = A >= critA;
B = zeros(size(A));
patches =zeros(1,ceil(numel(A)/10));
ipatch = 1;
[siz1,siz2] = size(A);
newndx=zeros(8,2);
ndx1=zeros(100,2);
%nx1=1;
for i = 1:size(A, 1)
   for j = 1:size(A, 2)
      if C(i, j) && B(i,j)==0
         B(i, j) = ipatch;
         nx1=1;
         ndx1(nx1,:) = [i,j];
         nx1=nx1+1;
         nx=1;
         newndx(nx,:)=[i,j];
         nx=nx+1;
         while nx>1
            ndx = newndx(1:nx-1,:);
            nx=1;
            for k = 1:size(ndx, 1)
               n1 = ndx(k, 1)-1; 
               n2 = ndx(k, 2);
               if (n1 > 0) && (n1 <= siz1) && (n2 > 0) && (n2 <= siz2) && C(n1, n2) && (B(n1, n2) == 0)
                  B(n1, n2) = ipatch;
                  newndx(nx,:) = [n1, n2]; 
                  nx=nx+1;
               end
               n1=n1+2;
               if  (n1 > 0) && (n1 <= siz1) && (n2 > 0) && (n2 <= siz2) && C(n1, n2) && (B(n1, n2) == 0)
                  B(n1, n2) = ipatch;
                  newndx(nx,:) = [n1, n2]; 
                  nx=nx+1;
               end
               n1 = ndx(k, 1); n2 = n2-1;
               if (n1 > 0) && (n1 <= siz1) && (n2 > 0) && (n2 <= siz2) && C(n1, n2) && (B(n1, n2) == 0)
                  B(n1, n2) = ipatch;
                  newndx(nx,:) = [n1, n2]; 
                  nx=nx+1;
               end
               n2=n2+2;
               if (n1 > 0) && (n1 <= siz1) && (n2 > 0) && (n2 <= siz2) && C(n1, n2) && (B(n1, n2) == 0)
                  B(n1, n2) = ipatch;
                  newndx(nx,:) = [n1, n2]; 
                  nx=nx+1;
               end
            end
            if nx>1
               ndx1(nx1:nx1+nx-2,:) =  newndx(1:nx-1,:);
               nx1=nx1+nx-1;
            end
         end
         s = nx1-1;
         %       B(sub2ind(siz, ndx1(:, 1), ndx1(:, 2))) = s;
         patches(ipatch) = s;
         ipatch = ipatch + 1;
      end
   end
end
patches = sort(patches(1:ipatch-1));
sizes = zeros(1,length(patches));
freq = zeros(1,length(patches));
i = 1;
j = 1;
while i <= length(patches)
   sizes(j) = patches(i);
   freq(j) = 1;
   i = i + 1;
   while i <= length(patches) && patches(i) == sizes(j)
      i = i + 1;
      freq(j) = freq(j) + 1;
   end
   j = j + 1;
end
if nargout == 0
   hfig=i_figure(100);
   plot(sizes, freq, 'o');
   try
      i_plotdefaults(hfig);
   catch
   end
   set(gca,'YScale','log');
   set(gca,'XScale','log');
   xlabel('size patch (cells)');
   ylabel('frequency');
else
   sizes1 = sizes;
   freq1 = freq;
   B1 = B;
end




