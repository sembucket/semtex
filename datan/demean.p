program demean(input, output);
{
    Remove mean values from numeric columnar input.

    Can be used as a filter.

    $Id$
}
const
    MAXCOL  = 32;

type
    rowofinput  = array [1..MAXCOL] of real;

    filehandle  = array [1..MAXCOL] of file of real;

var
    numcols, numrows,
    i, j            : integer;
    datum           : real;
    store           : filehandle;
    first           : rowofinput;
    meanvalue       : rowofinput;


function readfirst( {from}  var infile  : text;
                    {read}  var first   : rowofinput) : integer;
var
    ncols    : integer;
begin
    ncols := 1;
    while (not eoln(infile)) and (ncols <= MAXCOL) do begin
        read(infile, first[ncols]);
        ncols := ncols + 1
        end;
    ncols := ncols - 1;
    if (not eoln(infile)) then begin
        writeln(StdErr,
            'demean: more columns of input than available storage (',
            MAXCOL:1, ')');
        halt
        end;
    readln(infile);
    readfirst := ncols
end;



begin
    numcols := readfirst(input, first);
    numrows := 1;

    for i := 1 to numcols do begin
        rewrite(store[i]);
        meanvalue[i] := first[i];
        store[i]^ := first[i];
        put(store[i])
        end;
    
    while (not eof(input)) do begin
        for i := 1 to numcols do begin
            read(input, datum);
            meanvalue[i] := meanvalue[i] + datum;
            store[i]^ := datum;
            put(store[i])
            end;
        readln(input);
        numrows := numrows + 1
        end;
    
    for i := 1 to numcols do begin
        meanvalue[i] := meanvalue[i] / numrows;
        reset(store[i])
        end;
    
    for j := 1 to numrows do begin
        for i := 1 to numcols do begin
            datum := store[i]^;
            get(store[i]);
            datum := datum - meanvalue[i];
            write(output, ' ', datum:17)
            end;
        writeln(output)
        end
end.
