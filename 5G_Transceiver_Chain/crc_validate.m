function remainder=crc_validate(Rxed_Data,Gen)

value=Rxed_Data; 
array=value(1:length(Gen));

if array(1)==1
   remainder=double(xor(array,Gen));
   remainder=remainder(2:end);
else
    remainder=double(xor(array,zeros(1,length(Gen))));
    remainder=remainder(2:end);
end 
% through this if-else, first bit of data is checked, if it is 0, it is
% skipped, else it is xor with generator polynomial


for i=(length(Gen)+1):(length(value))
    array=[remainder,value(i)];
    if array(1)==1
        remainder=double(xor(array,Gen));
        remainder=remainder(2:end);
    else
        remainder=double(xor(array,zeros(1,length(Gen))));
        remainder=remainder(2:end);
    end
end
end
% if remainder array is all 0's, then there was error else some error
% occured