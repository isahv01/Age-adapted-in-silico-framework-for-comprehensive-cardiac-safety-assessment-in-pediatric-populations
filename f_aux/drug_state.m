%%
function [efecto_farmaco, namedrug] = drug_state(Num, d)
%Num = ID del fármaco a simular, 0 si control.
    if Num == 0
            efecto_farmaco = [1 1 1 1 1 1 1];
            namedrug = 'Control';
    else 
            %lectura de datos
            data = readtable('data_drugs2022_child.xlsx');
            %nombres_farmacos = data(:,1);      %nombres de los fármacos
            data_farmaco = data(Num,:);         %efecto de los fármacos
            
            %extracción de datos sobre farmaco a simular
            namedrug = char(data_farmaco.Drug_2022); %nombre
            data_farmaco = table2array(data_farmaco(1, 3:end)); %IC50 - h (7 corrientes)
            
            efecto_farmaco = drug_effect(data_farmaco, d);
    end
end
