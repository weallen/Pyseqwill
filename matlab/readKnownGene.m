function [kgdata, kgxrefdata] = readKnownGene( kgfname, kgxreffname )
%READ Summary of this function goes here
%   Detailed explanation goes here

    kg_file = fopen(kgfname);
    kgxref_file = fopen(kgxreffname);
    kgdata = textscan(kg_file, '%s %s %c %u8 %u8 %u8 %u8 %u8 %s %s %s %s');
    kgxrefdata = textscan(kgxref_file, '%s %s');
    fclose(kg_file);
    fclose(kgxref_file);    
end

