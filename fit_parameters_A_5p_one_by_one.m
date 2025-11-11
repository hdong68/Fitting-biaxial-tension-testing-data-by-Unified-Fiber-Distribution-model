
%%
clc; clear; 
clear all


Ns = 10;  % number of sample
directory = strings(Ns,1);

N_protocal=5;

%path_input='C:\Users\hdong68\Doc_Research\DocTBAD\1a_2022_TBAD_Tissue_Testing_Biaxial\B1_Chronic_Data';

%path_input='C:\Users\Arshiya\Documents\B1_Chronic_Data';

path_input='.\Selected_Data_5\';

% change 1 input

directory(1) = [path_input];
%xlsxfile(1) = [directory '\Sp4_cycle_7_m4_out - Copy.xls'];

directory(2) = [path_input '\2_adventitia_Feb_2_2022_LES020.041_FrT_Feb4_2022_Max_load430_S140kPa\analyzed\Selected_for_fitting_5p\'];
%xlsxfile(2) = [directory '\Sp4_cycle_3_m4_out - Copy.xls'];

directory(3) = [path_input '\3_adventitia_Feb_8_2022_LES115.247_FrT_Max_load470_S180kPa\analyzed\Selected_for_fitting_5p\'];
%xlsxfile(3) = [directory '\p4_cycle_4_m4_out - Copy.xls'];

directory(4) = [path_input '\4_adventitia_Feb_11_2022_1091_Max_load420_S120kPa\analyzed\Selected_for_fitting_5p\'];
%xlsxfile(4) = [directory '\Sp12_cycle_7_m12_out - Copy.xls'];

directory(5) = [path_input '\5_adventitia_Feb_15_2022_LES126.267_FrT_Max_load560_S210kPa\analyzed\Selected_for_fitting_5p\'];
%xlsxfile(5) = [directory '\Sp4_cycle_9_m4_out.xls'];

directory(6) = [path_input '\6_adventitia_Mar_18_2022_LES149.313_FrT_MaxLoad_300\analyzed\Selected_for_fitting_5p\'];
%xlsxfile(6) = [directory '\Sp12_cycle_5_m12_out - Copy.xls'];

directory(7) = [path_input '\7_11_adventitia_Apr5_2022_LES017.035_FrT_Outer_MaxLoad450\analyzed\Selected_for_fitting_5p\'];
%xlsxfile(7) = [directory '\Sp4_cycle_3_m4_out - Copy.xls'];

directory(8) = [path_input '\8_12_adventitia_Apr8_2022_LES033.081_FrT_MaxLoad420\analyzed\Selected_for_fitting_5p\'];
%xlsxfile = [directory '\Sp4_cycle_2_m4_out - Copy.xls'];

directory(9) = [path_input '\9_13_adventitia(guess)_Apr12_2022_LES007.015_FrT_MaxLoad320\analyzed\Selected_for_fitting_5p\'];
%xlsxfile = [directory '\Sp4_cycle_2_m4_out - Copy.xls'];

directory(10) = [path_input '\11_adventitia_15_Apr22_2022_LES118.253_FrT_Outer_MaxLoad450\analyzed\Selected_for_fitting_5p\'];
%xlsxfile = [directory '\Sp4_cycle_2_m4_out - Copy.xls'];




%%


ck_all=zeros(Ns,4); 
Rsq_all=zeros(Ns,6);

% change 2 initial guess

%ck0_all=xlsread('.\ck_Rsq_all_3p.xlsx','B2:E11');

%change 3 output figure 
path='.\Figures5P_one_by_one\';

%%

for jj=1:1%Ns %
    
    close all
    %jj=4;
    
    jj

   
    filePath=char(directory(jj));            
    
    filesAndFolders = dir(filePath);     % Returns all the files and folders in the directory

    filesInDir = filesAndFolders(~([filesAndFolders.isdir]));  % Returns only the files in the directory  
    numOfFiles = length(filesInDir);
    
    
    
    clear input;
    input=strings(N_protocal,1);

    for ii=1:N_protocal 
        %ii
        input(ii)=[filePath filesInDir(2*ii).name];
    end
    
    
    %%
    E11=[];E22=[];
    S11=[];S22=[];

    close all;

    % combine 7 protacols into one
    for ii=1:N_protocal     %2%

        data=xlsread(input(ii));
        % get the 2nd PK stress and Green Strain in x and y direction
        clear stress_X stress_Y strain_X strain_Y;   

       if jj==1 % || jj==5 || jj==9
            stress_X = data(:,13)*10^3;
            stress_Y = data(:,16)*10^3;
            strain_X= data(:,9);
            strain_Y= data(:,12);   
       else
            stress_X = data(:,16)*10^3;
            stress_Y = data(:,13)*10^3;
            strain_X= data(:,12);
            strain_Y= data(:,9);   
       end

        figure(1)
        plot(strain_X,stress_X,'b-');  
        hold on
        plot(strain_Y,stress_Y,'g-'); 

        xlabel('Strain');
        ylabel('Stress');

        E11=[E11;strain_X];
        E22=[E22;strain_Y];
        S11=[S11;stress_X];
        S22=[S22;stress_Y];

    end

    
    
    image_output = [path num2str(jj) 'A.tif'];
    print('-dtiff', image_output)
    % fit parameters for 7 protocols

    fun=@(ck)ComputeStreeError_UFD_A(ck,E11,S11,E22,S22);

    options = optimoptions('lsqnonlin','FunctionTolerance',1e-10);

    ck0=[50000,25000,10,0];
    lb=[0,0,0,0];
    ub=[1e10, 10^10, 1000, 1];

    [ck resnorm,residual]= lsqnonlin(fun,ck0,lb,ub,options);

    ck_kPa=[ck(1)/1000 ck(2)/1000 ck(3) ck(4)];
    
    ck_all(jj,:)=ck_kPa;

    %% compare fitting results

    clear S11_th S22_th


    lambda11=sqrt(1+2*E11);
    lambda22=sqrt(1+2*E22);

    I4=lambda11.^2;
    I4p=lambda22.^2; % perpendicular

    E_par=(1-ck(4))*(I4-1); % parallele
    E_per=ck(4)*(I4p-1); % perpendicular

    W1=ck(1)/2;
    W4=ck(2)*(1-ck(4))*(E_par).*exp(ck(3)*(E_par.^2+E_per.^2));
    W4p=ck(2)*(ck(4))*(E_per).*exp(ck(3)*(E_par.^2+E_per.^2));

    S11_th=ck(1)*(1-1./(lambda11.^4.*lambda22.^2))+2*W4;
    S22_th=ck(1)*(1-1./(lambda11.^2.*lambda22.^4))+2*W4p;



    figure(3)
    plot(E11,S11,'bo',E11,S11_th,'rx');  
    hold on

    xlabel('E11');
    ylabel('S11');
    
    %path='.\Figures7P_i3\';
    
    image_output = [path num2str(jj) 'B.tif'];
    print('-dtiff', image_output)

     figure(4)
    plot(E22,S22,'bo',E22,S22_th,'rx'); 

    xlabel('E22');
    ylabel('S22');

    image_output = [path num2str(jj) 'C.tif'];
    print('-dtiff', image_output)

    S11_MPa=S11/10^6;
    S22_MPa=S22/10^6;
    S11_th_MPa=S11_th/10^6;
    S22_th_MPa=S22_th/10^6;




    %CALCULATE R^2---------------------------------------------
    %load fit

    % S11

    SStot=0;
    SSres=0;
    for aa=1:length(S11)
        SStot=SStot+(S11(aa)-mean(S11)).^2;
        SSres=SSres+(S11(aa)-S11_th(aa)).^2;
    end
    Rsq11=1-(SSres/SStot);



    % S22
    SStot22=0;
    SSres22=0;
    for aa=1:length(S22)
        SStot22=SStot22+(S22(aa)-mean(S22)).^2;
        SSres22=SSres22+(S22(aa)-S22_th(aa)).^2;
    end
    Rsq22=1-(SSres22/SStot22);



    % R-squre of S11 and S22 together

    Rsq_Total=(Rsq11+Rsq22)/2.0;

    
    %CALCULATE root mean squre error RMSE---------------------------------------------

    % S11
    Res_Sum11=0;
    for aa=1:length(S11)
        Res_Sum11=Res_Sum11+(S11(aa)-S11_th(aa)).^2;
    end

    MSE11=Res_Sum11/length(S11);
    RMSE11=sqrt(MSE11)/1000;


    Res_Sum22=0;
    for aa=1:length(S22)
        Res_Sum22=Res_Sum22+(S22(aa)-S22_th(aa)).^2;
    end

    MSE22=Res_Sum22/length(S22);
    RMSE22=sqrt(MSE22)/1000;

    RMSE1122=(RMSE11+RMSE22)/2.0;

    Rsq_all(jj,:)=[Rsq11 Rsq22 Rsq_Total RMSE11 RMSE22 RMSE1122];

    %%
    E11S11=[E11 S11 S11_th];
    E22S22=[E22 S22 S22_th];


end

%%

%output for all sample

%change 4 output excel 
   
Unit={'c (kPa)', 'k1 (kPa)',	'k2 (-)',	'zeta (-)',	'Rsq11',	'Rsq22',	'R^2 (total)',	'RMSE11', 'RMSE22','RMSE1122'};

xlswrite('.\ck_Rsq_all_5p_i3p.xlsx', Unit,1,'B1');

 ck_Rsq_all=[ck_all Rsq_all];
ck_Rsq_all_c=num2cell(ck_Rsq_all);

xlswrite('.\ck_Rsq_all_5p_i3p.xlsx', ck_Rsq_all_c,1,'B2');




