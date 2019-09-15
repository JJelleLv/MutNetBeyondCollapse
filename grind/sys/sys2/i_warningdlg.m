function h=i_warningdlg(id,varargin)
w= warning('query',id);
if strcmp(w.state,'on')
   warning(id, varargin{:});
   [msg,msgid]=lastwarn;
   msg=['Warning: ' msg];
   h1=msgbox(striphtml(msg,1),['Warning ' msgid],'warn','non-modal');
   if nargout>0
       h=h1;
   end

end

