classdef varunit
   %VARUNIT class to implement analysis of units, you can evaluate expressions with
   %the units. Error messages will appear if the units are not matching.
   %SI units (and some other) can be converted to eachother.
   %
   %Usage:
   %   a=varunit('g/m3')
   %   b=varunit('yr');
   %   a/b returns
   %      g m-3 yr-1
   %   a^2/(b+1) returns
   %      g2 m-6 yr-1
   %
   % Note: short-cut and (&&) and short-cut or (||) cannot be used.
   % vector mode only limiited
   % Some tricks:
   % degreeC = °C
   % mum = µm
   %
   properties
      subunits = {}
      powers  = [];
      convfactor = 1;
   end
   
   properties (Dependent)
      unit
   end
   
   properties (Constant)
      prefixes = {1e24    'Y'
      1e21    'Z'
      1e18    'E'
      1e15    'P'
      1e12    'T'
      1e09    'G'
      1000000    'M'
      1000    'k'
      100    'h'
      10    'da'
      1    ''
      0.1000    'd'
      0.0100    'c'
      1e-03    'm'
      1e-06    char(181)
      1e-06    'mu'
      1e-09    'n'
      1e-12    'p'
      1e-15    'f'
      1e-18    'a'
      1e-21    'z'
      1e-24    'y'};
      
      base_units={'m';'kg';'s';'A';'K';'mol';'cd'};
      
      derived_units={'rad'      1         '-'  'radian'
      'sr'       1        '-'   'sterradian'
      'Hz'       1        's-1'  'Herz, frequency'
      'N'        1        'kg m s-2'  'Newton, force'
      'Pa'       1        'kg m-1 s-2'  'Pascal, pressure'
      'J'        1        'kg m2 s-2'  'Joules, energy'
      'W'        1        'kg m2 s-3'  'Watt, power'
      'C'        1        's A'     'Coulomb, amount of electicity'
      'V'        1        'kg m2 s-3 A-1'  'Volt, electric potential difference'
      'F'        1        'kg-1 m-2 s4 A2'  'Farad, capacitance'
      'Omega'    1        'kg m2 s-3 A-2'  'Ohm electric resistance'
      'S'        1        'kg-1 m-2 s3 A2' 'Siemens electrical conductance'
      'Wb'       1        'kg m2 s-2 A-1'  'Weber magnetic flux'
      'T'        1        'kg s-2 A-1'     'Tesla magnetic field strength'
      'H'        1        'kg m2 s-2 A-2'  'Henry, inductance'
      'lm'       1        'cd'         'lumen, luminous flux'
      'lux'      1        'cd m-2'  'lux,illuminance'
      'Bq'       1        's-1'    'becquerel, radioactivity'
      'Gy'       1        'm2 s-2'  'gray, absorbed dose'
      'Sv'       1        'm2 s-2'   'sievert, equivalent dose'
      'kat'      1        's-1 mol'  'katal, catalytic activity'
      };
      
      other_units = {
      char(176)      2*pi/360     'rad'  'degree'
      '#'        1        '-'  'number'
      'h'        3600         's'      'hour'
      'L'        1/1000        'm3'   'liter'
      'min'        60         's'     'minute'
      'hr'        3600         's'    'hour'
      'd'       24*3600         's'  'day'
      'day'       24*3600         's'  'day'
      'y'       365.25*24*3600     's'  'year'
      'yr'      365.25*24*3600     's'  'year'
      'inch'    0.0254       'm'   'inch'
      'ft'      0.3048       'm'    'feet'
      'yd'      0.9144       'm'   'yard'
      'mi'     1609.344      'm'  'mile'
      'a'       365.25*24*3600     's' 'year'
      'l'        1/1000        'm3'  'liter'
      'bar'      1000*100000 'g m-1 s-2' 'bar, pressure'
      'ha'        100*100        'm2'  'hectare'
      [char(176) 'C']       [1  273.15] 'K'      'degrees Celcius'
      [char(176) 'F']       [5/9 459.67] 'K'         'degrees Fahrenheit'
      'kWh'      3.6e9     'g m2 s-2' 'kilowatt hour'
      'M'        1000      'mol m-3'  'Molair, concentration'
      'cal'      4186.8   'g m2 s-2' 'calorie'
      'eV'   1.60217656535E-19 'kg m2 s-2' 'electronVolt'
      };
   end
   
   
   
   methods
      
      function obj = varunit(s, s1)
         if nargin == 2
            obj.subunits = s;
            obj.powers = s1;
         elseif nargin==1&&iscell(s)
            for i = length(s):-1:1
               obj(i) = varunit(s{i});
            end
         elseif nargin==1&&isa(s, 'varunit')
            obj.subunits = s.subunits;
            obj.powers = s.powers;
            obj.convfactor = s.convfactor;
         elseif nargin==1&&ischar(s)
            obj.unit = s;
         end
      end
      function obj = set.unit(obj, s)
         s = strtrim(s);
         obj.subunits = {};
         obj.powers = [];
         obj.convfactor = 1;
         if ~isempty(s)&&~strcmp(s, '[]')
            f = strfind(s, '[');
            if ~isempty(f)
               f2 = strfind(s, ']');
               for i = length(f):-1:1
                  s = s([1:f(i) - 1, f2(i) + 1:end]);
               end
            end
            if strcmp(s, '-')
               obj.powers = -99;
               return;
            end
            if length(s) > 6
               s = strrep(s, 'degree', char(176));
               s = strrep(s, 'micro', char(181));
            end
            %        s=strrep(s,'mu',char(181));
            f = regexp(s, '\<mu');
            for i = length(f):-1:1
               s = [s(1:f(i) - 1) char(181) s(f(i) + 2:end)];
            end
            s=strrep(s,'oC',[char(176) 'C']);
            s=strrep(s,'oF',[char(176) 'F']);
            lett=isletter(s)|s=='#'|s==char(176)|s=='µ';
            startw = find(diff([0 lett]) > 0);
            if isempty(startw)
               return;
            end
            endw = find(diff([lett 0]) < 0);
            obj.subunits = cell(1, length(startw));
            obj.powers = ones(1, length(startw));
            hasdiv= any(s(1:startw(1) - 1) == '/');
            hasdiv1 = false;
            hasbrack=any(s(1:startw(1) - 1) == '(');
            p = s(1:startw(1) - 1);
            p1=p(p~=' '&p~='.'&p~='*'&p~='/'&p~='^'&p~='('&p~=')');
            fac = str2num(p1);  %#ok<ST2NM>
            if ~any(isnan(fac))&&~any(isempty(fac))
               obj.convfactor = fac;
            end
            for i = 1:length(startw)
               obj.subunits{i} = s(startw(i):endw(i));
               if i < length(startw)
                  p = s(endw(i) + 1:startw(i + 1) - 1);
               else
                  p = s(endw(i) + 1:end);
               end
               p1=p(p~=' '&p~='.'&p~='*'&p~='/'&p~='^'&p~='('&p~=')');
               if isempty(p1)
                  apow = 1;
               else
                  apow = str2num(p1);  %#ok<ST2NM>
                  if all(isnan(apow))&&all(isempty(apow))
                     apow = 1;
                     warning('varunit:power','power neglected as it was no valid number')
                  end
               end
               if (hasdiv&&~hasbrack)||(hasbrack&&xor(hasdiv, hasdiv1))
                  apow = -apow;
               end
               obj.powers(i) = apow;
               hasbrack=hasbrack&&~any(p==')');
               if ~hasbrack
                  hasdiv=any(p == '/');
               else
                  hasdiv1=any(p == '/');
               end
               hasbrack=hasbrack||any(p=='(');
            end
         end
      end
      function s = get.unit(obj)
         if isundefined(obj)
            s = '[]';
         elseif isunitless(obj)
            s = '-';
         else
            s = sprintf('%s^%%g ', obj.subunits{:});
            s = sprintf(s, obj.powers);
            s=strrep(s,'^1 ',' ');
            s=strrep(s,'^','');
         end
         s = strtrim(s);
         if length(obj.convfactor) > 1||abs(obj.convfactor-1) > 1E-6
            if length(obj.convfactor) == 2
               s = sprintf('%g (+%g) %s', obj.convfactor, s);
            else
               s = sprintf('%g %s', obj.convfactor, s);
            end
            % s=sprintf('%g %s',obj.convfactor,s);
         end
      end
      function obj1 = simplify(obj)
         %simplify the unit
         if isunitless(obj)
            obj1 = obj;
            return;
         end
         usubunits = unique(obj.subunits);
         if length(usubunits) < length(obj.subunits)
            obj1 = varunit(obj);
            for i = 1:length(usubunits)
               ndx = find(strcmp(obj1.subunits, usubunits{i}));
               if length(ndx) > 1
                  apow = sum(obj1.powers(ndx));
                  if apow ~= 0
                     obj1.powers(ndx(1)) = apow;
                     ndx = ndx(2:end);
                  end
                  obj1.powers(ndx) = [];
                  obj1.subunits(ndx) = [];
               end
            end
         else
            obj1 = obj;
         end
         %unwritten rules for sorting subunits
         %
         units={'kg','g','s','d','mg','L','m','mol','µg','l','ha'};
         order = [1     1   4   4   1    2   2   1     1    2 , 2   ];
         order1 = zeros(size(obj1.subunits));
         for i = 1:length(obj1.subunits)
            ndx = strcmp(obj1.subunits{i}, units);
            if ~any(ndx)
               order1(i) = 10;
            else
               order1(i) =  order(ndx);
            end
         end
         [~, ndx] = sort(order1);
         obj1.subunits = obj1.subunits(ndx);
         obj1.powers = obj1.powers(ndx);
      end
      
      
      function [res, changed] = solveoper(objs, oper, constants, eqns,pos)
         function errormsg(varargin)
            if ~isempty(eqns)
               s = eqns.errorat(pos);
            end
            error('varunit:solveroper','%s',s,sprintf(varargin{1},varargin{2:end}));
         end
         function  fixedunitfcn(units)
            for k = 1:length(res)
               if  k < length(units)&&~isempty(units{k})
                  if notdefined(k)
                     res(k).unit  = units{k};
                     changed(k) = true;
                  else
                     res1 = varunit(units{k});
                     if ~res1.sameunit(res(k))
                        errormsg('Argument of "%s" has wrong unit (is [%s] should be [%s])', oper, res(k).unit, res1.unit);
                     end
                  end
               end
            end
         end
         function sameunitsfnc(mask, defunit)
            if length(res) < length(mask)
               mask = mask(1:length(res));
            end
            ndx1 = notdefined;
            ndx1(~mask) = true;
            f1 = find(~ndx1, 1);
            ndx1(~mask) = false;
            if ~isempty(f1)
               res(ndx1) = res(f1);
               changed(ndx1) = true;
            end
            ndx1 = notdefined & ~mask;
            if nargin==2&&any(ndx1)
               f = find(ndx1);
               for l = 1:length(f)
                  res(f(l)).unit = defunit;
                  changed(f(l)) = true;
               end
            end
            %check
            res1 = [];
            for l = 1:length(res)
               if mask(l)
                  if isempty(res1)
                     if ~res(l).isundefined
                        res1 = res(l);
                     end
                  elseif (~res(l).isundefined)&& ~res1.sameunit(res(l))
                     errormsg('Arguments of "%s" have different units ([%s] and [%s])', oper, res1.unit, res(l).unit)
                  end
               end
            end
         end
         if nargin < 4
            eqns = [];
         end
         res = objs;
         changed = false(size(objs));
         notdefined = objs.isundefined;
         isdefined = ~notdefined;
         isconst = ~isnan(constants);
         %         isinteger=round(constants)==constants;
         %          if ~any(isdefined)&&any(isinteger)
         %              for i=1:length(res)
         %                 if isinteger(i)&&~isdefined(i)
         %                      res(i).unit = '-';
         %                      changed(i) = true;
         %                 end
         %              end
         %              notdefined = res.isundefined;
         %              isdefined=~notdefined;
         %          end
         %   if any(notdefined)
         switch oper
          case {'+','-','''','.''','','ceil','floor','round','sum','mean','min',...
               'max','transpose','abs','real'} %including unary + and - all units the same
            sameunitsfnc([true true true]);
          case {'Dx','Dxx','Dy','Dyy'}
            sameunitsfnc([true true], '-')
          case 'iif'
            sameunitsfnc([true false true true], '-')
          case 'rednoise'
            %Tn=rednoise[t,T0,lambda,beta,iset,deltat)
            sameunitsfnc([true false true false true false false])
            if notdefined(4)
               res(4).unit = '-';
               changed(4) = true;
            elseif ~isunitless(res(4))
               if ~res(4).sameunit(res(1))
                  errormsg('Fourth arguments of "rednoise" is wrong ([%s] should be [-])', res(4).unit);
               end
            end
          case {'dwiener'} %result has unit of first arg
            sameunitsfnc([true true false false])
          case {'leftcells','rightcells','upcells','downcells', 'neighborcells','lag','repmat'} %result has unit of first arg
            sameunitsfnc([true true false false], '-')
          case {'>','>=','<','<=','==','~='} %result unitless
            sameunitsfnc([false true true], '-')
          case {'^','.^','sqrt'} %power note that a value is required
            if strcmp(oper, 'sqrt')
               constants(3) = 0.5;
            end
            if length(constants) < 3
               constants(3) = NaN;
            end
            if length(notdefined) > 2&&notdefined(3)
               res(3).unit = '-';
               notdefined(3) = false;
               changed(3) = true;
            end
            f = find(notdefined);
            f1 = find(isdefined);
            if length(f)==1&&isunitless(res(f1(1)))
               if notdefined(1)
                  res(1).unit = '-';
                  changed(1) = true;
               end
               if notdefined(2)
                  res(2).unit = '-';
                  changed(2) = true;
               end
            elseif length(f)==1&&f==2&&~isnan(constants(3))
               res(2) = res(1)^(1 / constants(3));
               changed(2) = true;
            elseif length(f)==1&&f==1&&~isnan(constants(3))&&(constants(3)~=0)
               res(1) = res(2)^constants(3);
               changed(1) = true;
            elseif isempty(f)
               if ~isnan(constants(3))
                  r = res(2)^constants(3);
                  if ~r.sameunit(res(1))
                     errormsg('Units of arguments of "^" are wrong ([%s]~=[%s]^%g)', res(1).unit, res(2).unit, constants(3));
                  end
               elseif ~isunitless(res(3))
                  errormsg('Power should be dimensionless (is now[%s])', res(3).unit);
               end
            end
          case {'sin','cos','tan'} %argument in radians
            fixedunitfcn({'-','rad'});
          case {'asin','acos','atan'} %res in radians
            fixedunitfcn({'rad','-'});
          case {'rand','randn','randlogn'}
            fixedunitfcn({'','-','-'});
          case {'exp','log','ln','log10','&','&&','|','||','~','not'} %unitless result and arguments
            fixedunitfcn({'-','-','-','-','-'});
          case {'*','.*'}  %multiply/divide
            f = find(notdefined);
            if isempty(f)
               r = res(2) * res(3);
               if ~r.sameunit(res(1))
                  errormsg('Units of arguments of "*" are wrong ([%s]~=[%s]*[%s])', res(1).unit, res(2).unit, res(3).unit);
               end
            elseif length(f) == 1
               switch f
                case 1
                  res(1) = res(2) * res(3);
                  changed(1) = true;
                case 2
                  res(2) = res(1) / res(3);
                  changed(2) = true;
                case 3
                  res(3) = res(1) / res(2);
                  changed(3) = true;
                otherwise
                  errormsg('Too many arguments for multiplication');
               end
            end
          case {'/','./'} %multiply/divide
            f = find(notdefined);
            if isempty(f)
               r = res(2) /  res(3);
               if ~r.sameunit(res(1))
                  errormsg('Units of arguments of "/" are wrong ([%s]~=[%s]/[%s])', res(1).unit, res(2).unit, res(3).unit);
               end
            elseif length(f) == 1
               switch f
                case 1
                  res(1) = res(2) / res(3);
                  changed(1) = true;
                case 2
                  res(2) = res(1) * res(3);
                  changed(2) = true;
                case 3
                  res(3) = res(2) / res(1);
                  changed(3) = true;
                otherwise
                  errormsg('Too many arguments for division');
               end
            end
          case {'size','length','numel','isnan','isa','ischar'}
            fixedunitfcn({'-'});
          case {'solver','mod','val','zeros','ones','eye','roundconv2',...
               'conv2','heaviside','implicitdisperse','prod','boxcartrain',...
               'boxcarinflow','parlookup','colon','horzcat','vertcat'}
            %skip suppress warning
          otherwise
            warning('varunit:solveroper','Operator/function "%s" is unknown',oper);
         end
         %      end
         changed(isconst) = false;
      end
      function res = convert(obj, d, a)
         if nargin == 2
            a = d;
            d = 1;
         end
         if ischar(a)
            a = varunit(a);
         end
         s2 = obj.basic;
         s1 = a.basic;
         if sameunit(s1, s2, true)
            %[K] = ([°F] + 459.67) × 5?/9	[°F] = [K] × 9/5 ? 459.67
            res1 = d;
            if length(s1.convfactor) > 1
               res1 = (res1 + s1.convfactor(2));
            end
            res1 = res1 * s1.convfactor(1) / s2.convfactor(1);
            if length(s2.convfactor) > 1
               res1 = (res1 - s2.convfactor(2));
            end
         else
            res1 = NaN;
         end
         if nargout == 0
            if isnan(res1)
               if ~s2.isbasic
                  fprintf('"%s" is unknown\n', char(obj));
               elseif ~s1.isbasic
                  fprintf('"%s" is unknown\n', char(a));
               else
                  fprintf('"%s" and "%s" do not match\n', char(a), char(obj));
               end
            else
               fprintf('%s %s = %s %s\n', mat2str(d), char(a), mat2str(res1), char(obj));
            end
         else
            res = res1;
         end
      end
      
      function res = isbasic(a)
         res = true;
         for i = 1:length(a.subunits)
            if ~any(strcmp(a.subunits{i}, a.base_units))
               res = false;
            end
         end
      end
      
      function res = isunitless(a)
         res=isempty(a.subunits)&&length(a.powers)==1&&a.powers==-99;
      end
      
      function res = isundefined(a)
         for i = length(a):-1:1
            res(i) = isempty(a(i).subunits)&&isempty(a(i).powers);
         end
      end
      
      function res = sameunit(a, b, ignore_convfactor)
         if nargin < 3
            ignore_convfactor = false;
         end
         if ischar(b)
            b = varunit(b);
         end
         res =  isunitless(a)&&isunitless(b);
         if ~res
            res=length(a.subunits)==length(b.subunits)&&all(a.powers==b.powers)&&length(union(a.subunits, b.subunits))==length(a.subunits)&&all(strcmp(a.subunits, b.subunits))&& ...
               length(a.convfactor)==length(b.convfactor)&&all(a.convfactor==b.convfactor);
            if ~res
               %check for synonyms
               a1 = a.basic;
               b1 = b.basic;
               if ignore_convfactor
                  a1.convfactor = 1;
                  b1.convfactor = 1;
               end
               if length(a1.subunits) ~= length(b1.subunits)
                  res = false;
                  return;
               end
               [~, ndx1] = sort(a1.subunits);
               [~, ndx2] = sort(b1.subunits);
               a1.subunits = a1.subunits(ndx1);
               b1.subunits = b1.subunits(ndx2);
               a1.powers = a1.powers(ndx1);
               b1.powers = b1.powers(ndx2);
               res=length(a1.subunits)==length(b1.subunits)&&all(a1.powers==b1.powers)&&length(union(a1.subunits, b1.subunits))==length(a1.subunits)&&all(strcmp(a1.subunits, b1.subunits))&& ...
                  length(a1.convfactor)==length(b1.convfactor)&&all(a1.convfactor==b1.convfactor);
            end
         end
      end
      
      %operator overloads
      function obj = plus(a, b)
         obj = plusoperator('+', a, b);
      end
      function obj = minus(a, b)
         obj = plusoperator('-', a, b);
      end
      function obj = lt(a, b)
         obj = gtoperator('<', a, b);
      end
      function obj = le(a, b)
         obj=gtoperator('<=', a, b);
      end
      function obj = gt(a, b)
         obj = gtoperator('>', a, b);
      end
      function obj = ge(a, b)
         obj=gtoperator('>=', a, b);
      end
      function obj = ne(a, b)
         obj=gtoperator('>=', a, b);
      end
      function obj = eq(a, b)
         obj=gtoperator('==', a, b);
      end
      function obj = and(a, b)
         %and  & operation for units
         obj = nondimfunction('&', a, b);
      end
      function obj = or(a, b)
         %or | operation for units
         obj = nondimfunction('|', a, b);
      end
      function obj = not(a)
         obj = nondimfunction('~', a);
      end
      function obj = uplus(a)
         obj = varunit(a);
      end
      function obj = uminus(a)
         obj = varunit(a);
      end
      function obj = times(a, b)
         obj = timesoperator(false, a, b);
      end
      function obj = mpower(a, b)
         obj = poweroperator(a, b);
      end
      function obj = sqrt(a)
         obj = poweroperator(a, 0.5);
      end
      function obj = power(a, b)
         obj = poweroperator(a, b);
      end
      function obj = mtimes(a, b)
         obj = timesoperator(false, a, b);
      end
      function obj = exp(a)
         obj = nondimfunction('exp(x)', a);
      end
      function obj = sin(a)
         obj = nondimfunction('sin(x)', a);
      end
      function obj = cos(a)
         obj = nondimfunction('cos(x)', a);
      end
      function obj = tan(a)
         obj = nondimfunction('tan(x)', a);
      end
      function obj = transpose(a)
         obj = a;
      end
      function obj = ctranspose(a)
         obj = a;
      end
      function obj = colon(varargin)
         if isnumeric(varargin{1})
            obj = varunit('');
         else
            obj = varargin{1};
         end
      end
      function obj = mrdivide(a, b)
         obj = timesoperator(true, a, b);
      end
      function obj = rdivide(a, b)
         obj = timesoperator(true, a, b);
      end
      function obj = ldivide(a, b)
         obj = timesoperator(true, b, a);
      end
      function obj = horzcat(a, b)
         obj = plusoperator('horzcat', a, b);
      end
      function obj = vertcat(a, b)
         obj = plusoperator('vertcat', a, b);
      end
      function obj = mldivide(a, b)
         obj = timesoperator(true, b, a);
      end
      function s = char(obj)
         if length(obj) > 1
            s = '';
            for i = 1:length(obj)
               s = sprintf('%s%s\n', s, obj(i).unit);
            end
            return;
         end
         s = obj.unit;
      end
      function disp(obj)
         disp(char(obj));
      end
      
      function [obj, unknowns] = basic(a)
         if length(a) > 1
            for i = length(a):-1:1
               obj(i) = a(i).basic;
            end
            return;
         end
         unknowns = {};
         obj = varunit(a);
         objs = cell(size(a.subunits));
         for i = 1:length(a.subunits)
            n = a.subunits{i};
            objs{i} = [];
            prefix = '';
            factor1 = 1;
            if length(n)>1&&~any(strcmp(n,{'cd','ha','Gy','Pa','yd','day', 'days'})) %some symbols have different interpretations
               %eg. ha hecto annellum or hectare (the latter is preferred)
               %cd = centiday makes no sense
               fprefix = find(strncmp(n, a.prefixes(:, 2), 1));
               if length(fprefix) > 1
                  longpref={'da','mu'};
                  fprefix2 = strncmp(n, longpref, 2);
                  if any(fprefix2)&&length(n) > 2
                     fprefix = find(strcmp(longpref{fprefix2}, a.prefixes(:, 2)));
                  else
                     fprefix = find(strcmp(n(1), a.prefixes(:, 2)));
                  end
               end
               if ~isempty(fprefix)
                  prefix = varunit.prefixes{fprefix, 2};
                  factor1 = varunit.prefixes{fprefix, 1};
               end
            end
            rest = n(length(prefix) + 1:end);
            if strcmp(rest, 'g')
               %for some reason kg is preferred by SI over g
               rest = 'kg';
               factor1 = factor1 / 1000;
            end
            if any(strcmp(rest, a.base_units))
               obj.subunits{i} = rest;
               if obj.powers(i)~=1&&factor1~=1
                  factor1 = factor1^obj.powers(i);
               end
               obj.convfactor = obj.convfactor * factor1;
            else
               f = find(strcmp(rest, a.derived_units(:, 1)));
               if isempty(f)
                  f = find(strcmp(n, a.derived_units(:, 1)));
                  if ~isempty(f)
                     factor1 = 1;
                  end
               end
               if ~isempty(f)
                  obj.subunits{i} = '';
                  objs{i} = varunit(a.derived_units{f, 3});
                  objs{i}.convfactor = factor1 * a.derived_units{f, 2};
                  if obj.powers(i) ~= 1
                     objs{i} = objs{i}^obj.powers(i);
                  end
               else
                  %last in preference are the other units
                  f = find(strcmp(rest, a.other_units(:, 1)));
                  if isempty(f)
                     f = find(strcmp(n, a.other_units(:, 1)));
                     if ~isempty(f)
                        factor1 = 1;
                     else
                        unknowns{end + 1} = n;
                     end
                  end
                  if ~isempty(f)
                     obj.subunits{i} = '';
                     objs{i} = varunit(a.other_units{f, 3});
                     objs{i}.convfactor = factor1 * a.other_units{f, 2};
                     if obj.powers(i) ~= 1
                        objs{i} = objs{i}^obj.powers(i);
                     end
                  end
               end
            end
         end
         ndx = ~strcmp(obj.subunits, '');
         obj.subunits = obj.subunits(ndx);
         obj.powers = obj.powers(ndx);
         for i = 1:length(objs)
            if ~isempty(objs{i})
               obj = obj * objs{i};
            end
         end
         obj = simplify(obj);
         if isempty(obj.powers)
            obj.powers(1) = -99;
         end
      end
   end
   methods (Access = private)
      
      function obj = nondimfunction(fun, varargin)
         for i = 1:length(varargin)
            if ~islogical(varargin{i})&&~isnumeric(varargin{i})&&~isunitless(varargin{i})
               error('varunit:plus','The arguments of the function "%s" should be dimensionless',fun);
            end
         end
         obj = varunit('-');
      end
      function obj = poweroperator(a, b)
         if isnumeric(a)
            obj = varunit;
         elseif isnumeric(b)
            obj = varunit(a);
            if ~isunitless(a)
               obj.powers = obj.powers * b;
            end
            obj.convfactor = obj.convfactor.^b;
         else
            if ~isunitless(b)&&~isundefined(b)
               error('varunit:power','Power must be unitless');
            end
            obj = varunit(a);
            warning('varunit:power','Cannot calculate variable power');
         end
         if ~isunitless(a)
            ndx=obj.powers ~= 0;
            obj.subunits = obj.subunits(ndx);
            obj.powers = obj.powers(ndx);
            if isempty(obj.powers)
               obj.powers = -99;
            end
         end
      end
      function obj = timesoperator(isdivision, a, b)
         if isnumeric(a)
            obj = varunit(b);
            % obj.powers = powfun(0, obj.powers); %dont understand this line
            return;
         end
         if isnumeric(b)
            obj = varunit(a);
            return;
         end
         if isunitless(a)&&isunitless(b)
            obj = a;
            return;
         end
         a = simplify(a);
         b = simplify(b);
         obj = varunit;
         %unsorted union
         obj.subunits = a.subunits;
         for i = 1:length(b.subunits)
            if ~any(strcmp(b.subunits{i}, a.subunits))
               obj.subunits(end + 1) = b.subunits(i);
            end
         end
         obj.powers = zeros(size(obj.subunits));
         for i = 1:length(obj.subunits)
            sub1 = strcmp(obj.subunits{i}, a.subunits);
            if any(sub1)
               pow1 = a.powers(sub1);
            else
               pow1 = 0;
            end
            sub2 = strcmp(obj.subunits{i}, b.subunits);
            if any(sub2)
               pow2 = b.powers(sub2);
            else
               pow2 = 0;
            end
            if isdivision
               obj.powers(i) = pow1(1) - pow2(1);
            else
               obj.powers(i) = pow1(1) + pow2(1);
            end
         end
         ndx=obj.powers ~= 0;
         obj.subunits = obj.subunits(ndx);
         obj.powers = obj.powers(ndx);
         if isdivision
            obj.convfactor = a.convfactor / b.convfactor;
         else
            obj.convfactor = a.convfactor * b.convfactor;
         end
         if isempty(obj.powers)
            obj.powers = -99;
         end
         obj = obj.simplify;
      end
      function obj = gtoperator(oper, a, b)
         obj = varunit('-');
         if ~(isnumeric(a)||isnumeric(b))
            if ~(isundefined(a)||isundefined(b))&&~sameunit(a, b)
               error('varunit:plus','The units should be the same for this operation "%s"',oper);
            end
         end
      end
      function obj = plusoperator(oper, a, b)
         if isnumeric(a)
            obj = varunit(b);
            return;
         end
         if isnumeric(b)
            obj = varunit(a);
            return;
         end
         if ~(isundefined(a)||isundefined(b))&&~sameunit(a, b)
            error('varunit:plus','The units should be the same for this operation "%s"',oper);
         end
         obj = varunit(a);
      end
   end
end




