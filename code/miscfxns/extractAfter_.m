function newStr = extractAfter_(str,startStrPos)
if( isnumeric(startStrPos) )
    if( startStrPos > numel(str) )
        error('Numeric position exceeds the number of characters.');
    else
        newStr = str(startStrPos+1:end);
    end
else
    s = strfind(str,startStrPos);
    if( isempty(s) )
        error('Conversion from <missing> to character vector is not supported.');
    else
        newStr = str(s(1)+numel(startStrPos):end);
    end
end
end