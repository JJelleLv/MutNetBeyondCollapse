% ---------------------------------------------------------------------------
function i_slider_Callback(obj,~,hand)
if nargin<3
  hfig=getparentfig(obj);
  hand= get(hfig,'userdata');
end

aval = round(get(hand.hSlid,'Value'));
old_val = get(hand.hSlid,'UserData');
ds = aval - old_val;

if (ds < 0)                                         % Slider moved down
    n = hand.NumRows - aval + 1;    d_col = hand.NumRows - aval;
    if (n+hand.MAX_ROWS-1 > hand.NumRows)             % Case we jumped into the midle zone
        adj = (n+hand.MAX_ROWS-1 - hand.NumRows);
        n = n - adj;    d_col = d_col - adj;
    end
    for i = n:min(n+hand.MAX_ROWS-1,hand.NumRows)   % Update positions
        for j = 1:hand.NumCol
            pos = hand.Edits_pos_orig{i,j};
            set(hand.hEdits(i,j),'pos',[pos(1) hand.arr_pos_y(i-d_col) pos(3:4)],'Visible','on')
        end
        if (~isempty(hand.checks))                  % If we have checkboxes
            pos = hand.Checks_pos_orig(i,:);
            set(hand.hChecks(i),'pos',[pos(1) hand.arr_pos_y(i-d_col)+3 pos(3:4)],'Visible','on')            
        end
    end
    if (i == get(hand.hEdits(hand.NumRows,1),'UserData')) % Bottom reached. Jump to there
        aval = 1;    set(hand.hSlid,'Value',aval)         % This also avoids useless UIs repositioning
    end
elseif (ds > 0)                                     % Slider moved up
    n = hand.NumRows - aval + 1;    k = hand.MAX_ROWS;
    if (n < hand.MAX_ROWS)                          % Case we jumped into the midle zone
        adj = (hand.MAX_ROWS - n - 0);
        n = n + adj;
    end
    for i = n:-1:max(n-hand.MAX_ROWS+1,1)         % Update positions
        for j = 1:hand.NumCol
            pos = hand.Edits_pos_orig{i,j};
            set(hand.hEdits(i,j),'pos',[pos(1) hand.arr_pos_y(k) pos(3:4)],'Visible','on')        
        end
        if (~isempty(hand.checks))                  % If we have checkboxes
            pos = hand.Checks_pos_orig(i,:);
            set(hand.hChecks(i),'pos',[pos(1) hand.arr_pos_y(k)+3 pos(3:4)],'Visible','on')            
        end
        k = k - 1;
    end
    set(hand.hEdits(n+1:end,1:end),'Visible','off')
    if (~isempty(hand.checks)),     set(hand.hChecks(n+1:end),'Visible','off');    end
    if (i == get(hand.hEdits(1,1),'UserData'))      % Reached Top. Jump to there
        set(hand.hSlid,'Value',hand.NumRows)          % This also avoids useless UIs repositioning
        aval = hand.NumRows;
    end
end
set(hand.hSlid,'UserData',aval)                      % Save old 'Value'
