function [hexOut] = dec2hexTwoComp(decIn)
%DEC2HEXTWOCOMP Creates the Hexadecimal two's complement from an input
%int16, which could be type cast as a double.
%   Input is assumed to be a int16, i.e., [-2^15, 2^15 - 1] = [-32768, 32767]

% Initialize the output. Ndat long and 4 characters wide for 16bit hex
Ndat = length(decIn);
hexOut = repmat(char(80), Ndat,4); % char(80) is arbitrary, just want a string

for idx = 1:Ndat
    if decIn(idx) >= 0
        % If number isn't negative, just get normal representation
        hexOut(idx,:) = dec2hex(decIn(idx), 4);
    else
        % If number is negative, get the 2's complement representation
        % Convert the absolute value to binary array
        binIn  = de2bi(abs(decIn(idx)), 16);

        % Flip each bit
        binOut = heaviside(.5 - binIn);

        % Turn back into dec and add one
        decOut = bi2de(binOut) + 1;

        % Turn into hex
        hexOut(idx,:) = dec2hex(decOut, 4);
    end
end
    
end