function efecto_farmaco = drug_effect(data_farmaco, d)
%data_farmaco = vector 1x15 (quitamos el nombre y la clase antes)
%d = valor que indica el nº de veces de concentracion de EFTPC

    %extraccion de datos por canal
    D = data_farmaco(1); %EFTPC
    %Ikr
    IC50 = data_farmaco(2);
    h = data_farmaco(3);
    factor_Ikr = (1+(d*D/IC50)^h)^-1;
    %INa
    IC50 = data_farmaco(4);
    h = data_farmaco(5);
    factor_INa = (1+(d*D/IC50)^h)^-1;
    %INaL
    IC50 = data_farmaco(6);
    h = data_farmaco(7);
    factor_INaL = (1+(d*D/IC50)^h)^-1;
    %ICaL
    IC50 = data_farmaco(8);
    h = data_farmaco(9);
    factor_ICaL = (1+(d*D/IC50)^h)^-1;
    %IKs
    IC50 = data_farmaco(10);
    h = data_farmaco(11);
    factor_Iks = (1+(d*D/IC50)^h)^-1;
    %IK1
    IC50 = data_farmaco(12);
    h = data_farmaco(13);
    factor_IK1 = (1+(d*D/IC50)^h)^-1;
    %Ito
    IC50 = data_farmaco(14);
    h = data_farmaco(15);
    factor_Ito = (1+(d*D/IC50)^h)^-1;

    %formar vector final
    efecto_farmaco = [factor_Ikr, factor_INa, factor_INaL, factor_ICaL, factor_Iks, factor_IK1, factor_Ito];
    efecto_farmaco(isnan(efecto_farmaco)) = 1; %los NaN los sustituye por 1. 

end