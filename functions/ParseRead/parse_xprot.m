function mrprot = parse_xprot(rawarr)
% parse XProtocol structure stored in a char array
%  E. Auerbach, CMRR, 2013

mrprot = [];
tagName = cell(10);

% first, make sure it is a char array
if ~ischar(rawarr), error('parse_xprot: input is not a char array!'); end

% find <XProtocol>{
xstart = strfind(rawarr,'<XProtocol>');

if (xstart)
    workarr = extractBraceString(rawarr(xstart:end)); % <XProtocol> {}
    mrprot = parse_loop(mrprot, workarr, tagName, 0);
else
    error('parse_xprot() failed to find <XProtocol>!');
end

%--------------------------------------------------------------------------
function mrprot = parse_loop(mrprot, workarr, tagName, level)
% this function is called recursively to parse the xprotocol

[tagStr, tagType, workarr] = findNextTag(workarr);
while (~isempty(workarr))

    % for parammap, add the name as another level of the current array name
    % and spawn off another copy of this function to deal with them
    if (strncmpi(tagType, 'ParamMap', 8))
        if (~isempty(tagStr))
            level = level + 1;
            tagName{level} = make_safe_fieldname(tagStr);
        end
        stubArr = extractBraceString(workarr);
        mrprot = parse_loop(mrprot, stubArr, tagName, level);
        workarr = workarr(length(stubArr)+1:end);
        if (~isempty(tagStr))
            level = level - 1;
        end

    % for paramarray, special processing is required
    % for now, though, just skip them since the most important ones are
    % still in the text mrprot and it would be a lot more work to process
    % them here
    elseif (strncmpi(tagType, 'ParamArray', 10))
        if (~isempty(tagStr))
            level = level + 1;
            tagName{level} = make_safe_fieldname(tagStr);
        end
        %fprintf('found array %s\n',build_tag_name(tagName,level+1));
        stubArr = extractBraceString(workarr);
        workarr = workarr(length(stubArr)+1:end);

%         %skip the <default> tag, which always seems to be there
%         [tagStr, tagType, stubArr] = findNextTag(stubArr);
%         if (~strcmpi(tagStr,'Default'))
%             error('parse_xprot(): unknown ParamArray format!')
%         end
%         
%         % now map the fields
%         [tagStr, tagType, stubArr] = findNextTag(stubArr);
%         
%         mrarray = parse_loop([], stubArr, [], 0);
%         remArr = extractBraceString(stubArr);
%         stubArr = stubArr(length(remArr)+1:end);

        if (~isempty(tagStr))
            level = level - 1;
        end

    % these are values that we know how to process, so do so
    elseif ((strncmpi(tagType, 'ParamBool', 9))     || ...
            (strncmpi(tagType, 'ParamLong', 9))     || ...
            (strncmpi(tagType, 'ParamDouble', 11))  || ...
            (strncmpi(tagType, 'ParamString', 11)))
        tagName{level+1} = make_safe_fieldname(tagStr);
        stubStr = extractBraceString(workarr);
        %fprintf('%s = %s\n',tagName{level+1},stubStr);
        value = get_xprot_value(stubStr, tagType);
        fields = {tagName{1:level+1}, value};
        mrprot = setfield(mrprot, fields{:});
        workarr = workarr(length(stubStr)+1:end);

    % we don't care about the things below, but acknowledge that we know
    % about them
    elseif ((strncmpi(tagType, 'ParamChoice', 11))      || ...
            (strncmpi(tagType, 'Pipe', 4))              || ...
            (strncmpi(tagType, 'Connection', 10))       || ...
            (strncmpi(tagType, 'Event', 5))             || ...
            (strncmpi(tagType, 'ParamFunctor', 12))     || ...
            (strncmpi(tagType, 'Method', 6))            || ...
            (strncmpi(tagType, 'ProtocolComposer', 16)) || ...
            (strncmpi(tagType, 'Dependency', 10))       || ...
            (strncmpi(tagType, 'ParamCardLayout', 15))  || ...
            (strncmpi(tagType, 'EVACardLayout', 13)))
        stubStr = extractBraceString(workarr);
        workarr = workarr(length(stubStr)+1:end);

        % ignore these also
    elseif ((strncmpi(tagStr, 'Name', 4))             || ...
            (strncmpi(tagStr, 'ID', 2))               || ...
            (strncmpi(tagStr, 'Comment', 7))          || ...
            (strncmpi(tagStr, 'Label', 5))            || ...
            (strncmpi(tagStr, 'Visible', 7))          || ...
            (strncmpi(tagStr, 'Userversion', 11)))
        % skip to end of line for these
        lend = strfind(workarr, char(10));
        workarr = workarr(lend+1:end);
    elseif (strncmpi(tagStr, 'EVAStringTable', 14))
        % skip brace string for this one
        stubStr = extractBraceString(workarr);
        workarr = workarr(length(stubStr)+1:end);
        
    % something unknown has happened if we reach this point, so throw a warning
    else
        fprintf('parse_xprot(): WARNING: found unknown tag %s (%s)\n',tagStr,tagType);
    end

    [tagStr, tagType, workarr] = findNextTag(workarr);
end

%--------------------------------------------------------------------------
function value = get_xprot_value(stubStr, tagType)
% stubStr contains the values of interest, but also might include
% modifiers. we will just ignore the modifiers, which seem to always be
% terminated by CR.

[tagStr, newtagType, remStr] = findNextTag(stubStr);
while (~isempty(remStr)) % found a tag
    if (~isempty(newtagType)) % not a modifier tag???
        error('parse_xprot()::get_xprot_value(): ERROR: found unknown tag!')
    else
        % these are the modifier tags we know about
        if ((strncmpi(tagStr, 'Precision', 9))      || ...
            (strncmpi(tagStr, 'LimitRange', 10))    || ...
            (strncmpi(tagStr, 'MinSize', 7))        || ...
            (strncmpi(tagStr, 'MaxSize', 7))        || ...
            (strncmpi(tagStr, 'Limit', 5))          || ...
            (strncmpi(tagStr, 'Default', 7))        || ...
            (strncmpi(tagStr, 'InFile', 6))         || ...
            (strncmpi(tagStr, 'Context', 7))        || ...
            (strncmpi(tagStr, 'Dll', 3))            || ...
            (strncmpi(tagStr, 'Class', 5))          || ...
            (strncmpi(tagStr, 'Comment', 7))        || ...
            (strncmpi(tagStr, 'Label', 5))        || ...
            (strncmpi(tagStr, 'Tooltip', 7))        || ...
            (strncmpi(tagStr, 'Visible', 7))        || ...
            (strncmpi(tagStr, 'Unit', 4)))
            % acknowledge that we know about these tags
        else
            fprintf('parse_xprot(): WARNING: found unknown modifier %s\n',tagStr);
        end
        
        % remove the line containing the tag
        lstart = strfind(stubStr, ['<' tagStr '>']);
        lend = strfind(stubStr(lstart+1:end), char(10));
        if (lstart > 1)
            stubStr = [stubStr(1:lstart-1) stubStr(lend+lstart+1:end)];
        else
            stubStr = stubStr(lend+1:end);
        end
    end

    [tagStr, newtagType, remStr] = findNextTag(stubStr);
end

if (strncmpi(tagType, 'ParamString', 11))
    value = getQuotString(stubStr);
    ok = true;
else
    stubStr = strrep(stubStr, char(10), ' '); % remove newlines
    if (strncmpi(tagType, 'ParamBool', 9))
        stubStr = strrep(stubStr, '"true"', '1');
        [value,ok] = str2num(stubStr); %#ok<ST2NM>
        value = (value ~= 0);
    elseif (strncmpi(tagType, 'ParamLong', 9))
        [value,ok] = str2num(stubStr); %#ok<ST2NM>
    elseif (strncmpi(tagType, 'ParamDouble', 11))
        [value,ok] = str2num(stubStr); %#ok<ST2NM>
    end
end

if (~ok)
    fprintf('WARNING: get_xprot_value failed: %s\n', stubStr);
    disp(value);
end


%--------------------------------------------------------------------------

% function tagStr = build_tag_name(tagName, level)
% 
% tagStr = [];
% 
% if (level)
%     for x=1:level
%         if (x > 1), tagStr = strcat(tagStr, '.'); end
%         tagStr = strcat(tagStr, tagName{x});
%     end
% end

%--------------------------------------------------------------------------

function stvar = make_safe_fieldname(tagStr)
% this function checks potential fieldnames and makes sure they are valid
% for MATLAB syntax, e.g. 2DInterpolation -> x2DInterpolation (must begin
% with a letter)

tagStr = strtrim(tagStr);

if (isletter(tagStr(1)))
    stvar = tagStr;
else
    stvar = strcat('x', tagStr);
end

if strfind(stvar, ';'), stvar = strrep(stvar, ';', '_'); end
if strfind(stvar, '@'), stvar = strrep(stvar, '@', '_'); end % VD13

%--------------------------------------------------------------------------

function stvar = extractBraceString(text)
% extracts string from within curly braces, including nested braces

tlen = length(text);
tstart = strfind(text,'{');

stvar = [];
if (tstart)
    tnests = 1;
    for x=tstart+1:tlen
        if (text(x) == '{')
            tnests = tnests + 1;
        elseif (text(x) == '}')
            tnests = tnests - 1;
        end
        if (tnests == 0)
            stvar = text(tstart+1:x-1);
            break;
        end
    end
end

%--------------------------------------------------------------------------

function [tagStr, tagType, remStr] = findNextTag(inStr)
% returns <tag> name and the remainder of the string following the tag.
% for e.g. <Tag>, returns tagStr='Tag', tagType=''
% for e.g. <ParamLong."Tag">, returns tagStr='Tag', tagType='ParamLong'
% if no tag is found, returns null strings

tagStr = [];
tagType = [];
remStr = [];

startPos = strfind(inStr,'<'); % look for start of tag
if (startPos)
    endPos = strfind(inStr,'>'); % look for end of tag
    if (endPos)
        % found complete tag
        fullTag = inStr(startPos+1:endPos-1);
        
        % now check for name/type
        dotPos = strfind(fullTag,'."');
        if (dotPos)
            tagStr = getQuotString(fullTag(dotPos+1:end));
            tagType = fullTag(1:dotPos-1);
        else
            tagStr = fullTag;
        end
        
        % return remainder
        remStr = inStr(endPos+1:end);
    end
end

%--------------------------------------------------------------------------

function stvar = getQuotString(text)
% extracts string between double quotes, e.g. "string"
%  also works with double-double quotes, e.g. ""string""

idx = strfind(text,'"');

if ( (length(idx) == 4) && (idx(1)+1 == idx(2)) && (idx(3)+1 == idx(4)) ) % double-double quotes
    stvar = text(idx(2)+1:idx(3)-1);
elseif (length(idx) >= 2) % double quotes, or ??? just extract between first and last quotes
    stvar = text(idx(1)+1:idx(end)-1);
else % malformed?
    stvar = text;
end
