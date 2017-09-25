function blocks = findContiguousNZ(A)

wrap       = [0, A, 0] ;
temp       = diff( wrap ~= 0 ) ;
blockStart = find( temp == 1 ) + 1 ;
blockEnd   = find( temp == -1 ) ;
blocks     = arrayfun( @(bId) wrap(blockStart(bId):blockEnd(bId)), ...
    1:numel(blockStart), 'UniformOutput', false ) ;

end