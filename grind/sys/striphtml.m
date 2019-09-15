function s1 = striphtml(s, nohyper)
if nargin < 2
   nohyper = 0;
end
s=strrep(s,'==>','==&gt;');
f1 = strfind(s, '<');
f2 = strfind(s, '>');
f3 = [strfind(s, '<a') strfind(s, '</a')];
s1 = s;
if (~isempty(f1))
   if (length(f1) ~= length(f2))
       disp(s);
     error('striphtml:grind','unequal number of html tags, replace > with &gt; and < with &lt;')
   else
      for i = length(f1):-1:1
         if nohyper || ~any(f3==f1(i))
            if f1(i) == 1
               s1 = s1(f2(i) + 1:length(s1));
            elseif f2(i) == length(f1)
               s1 = s1(1:f1(i) - 1);
            else
               s1 = [s1(1:f1(i) - 1) s1(f2(i) + 1:length(s1))];
            end
            f1 = strfind(s1, '<');
            f2 = strfind(s1, '>');
         else
            f1 = f1(1:end - 1);
            f2 = f2(1:end - 1);
         end
      end
   end
end
s1 = symreplace('acute;', s1);
s1 = symreplace('grave;', s1);
s1 = symreplace('uml;', s1);
s1 = symreplace('ring;', s1);
s1 = symreplace('tilde;', s1);
s1 = symreplace('slash;', s1);
% add new special characters if used
s1=strrep(s1,'&copy;','(c)');
s1=strrep(s1,'&mu;',char(181));
s1=strrep(s1,'&nbsp;',' ');
s1=strrep(s1,'&lt;','<');
s1=strrep(s1,'&gt;','>');

function s1 = symreplace(symb, s)
f1 = length(strfind(s, symb));
if (f1 > 0) && (f1 == length(strfind(s, '&')))
   s1 = strrep(s, symb, '');
   s1=strrep(s1,'&','');
else
   s1 = s;
end
