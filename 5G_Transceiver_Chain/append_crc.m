function Tx_data=append_crc(message,polymial)

no_of_zeros=length(polymial)-1;
array=[message,zeros(1,no_of_zeros)]; %zeros are appended in the data to be txed

value=array(1:length(polymial));

if value(1)==1
   remainder=double(xor(value,polymial));
   remainder=remainder(2:end);
else
    remainder=double(xor(value,zeros(1,length(polymial))));
    remainder=remainder(2:end);
end 
% through this if-else, first bit of data is checked, if it is 0, it is
% skipped, else it is xor with generator polynomial

for i=(length(polymial)+1):(length(array))
    value=[remainder,array(i)];
    if value(1)==1
        remainder=double(xor(value,polymial));
        remainder=remainder(2:end);
    else
%         continue;
        remainder=double(xor(value,zeros(1,length(polymial))));
        remainder=remainder(2:end);

    end
end
crc=remainder;
Tx_data=[message,crc];  % CRC is appended in message block

end