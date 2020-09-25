%Plotter for results files of WENO Mott-Gurney solver
%Can plot from selection of single or multiple files, 
%just need to set single_file =true or false accordingly.

clear 

single_file = true;

[File,Path]=uigetfile('*.txt','MultiSelect','on');

N = numel(File);      %Counts # of files in cell array "File" that was outputted by uigetfile


if single_file
     
    format shortG   
    data= load ([Path File]);  

    x = data(:,1);
    p = data(:,2);
    E = data(:,3);
    Jp = data(:,4);

    h1 = plot(x,p);
    hold on
    title(File,'interpreter','latex','FontSize',16);
    xlabel('Position ($m$)','interpreter','latex','FontSize',14);
    ylabel({'Hole density ($1/m^3$)'},'interpreter','latex','FontSize',14);

    figure;
    h2 = plot(x,E);
    hold on
    title(File,'interpreter','latex','FontSize',16);
    xlabel('Position ($m$)','interpreter','latex','FontSize',14);
    ylabel({'Electric Field (V/m)'},'interpreter','latex','FontSize',14);


    figure;
    h3 = plot(x,Jp);
    hold on
    title(File,'interpreter','latex','FontSize',16);
    xlabel('Position ($m$)','interpreter','latex','FontSize',14);
    ylabel({'Current Density ($A/m^2$)'},'interpreter','latex','FontSize',14);

    hold off    
    
 else  %if only 1 file selected
    
     for  num=1:N         %Repeat loop for each file

        name= File(1,num);
        str=sprintf('%s', [Path name{1}]);          %makes str be the name of file (along with its path)

        format shortG   

        data= load (str);  

        x = data(:,1);
        p = data(:,2);
        E = data(:,3);
        Jp = data(:,4);

        h1 = plot(x,p);
        hold on
        title(name{1},'interpreter','latex','FontSize',16);
        xlabel('Position ($m$)','interpreter','latex','FontSize',14);
        ylabel({'Hole density ($1/m^3$)'},'interpreter','latex','FontSize',14);

        figure;
         h2 = plot(x,E);
        hold on
        title(name{1},'interpreter','latex','FontSize',16);
        xlabel('Position ($m$)','interpreter','latex','FontSize',14);
        ylabel({'Electric Field (V/m)'},'interpreter','latex','FontSize',14);


        figure;
        h3 = plot(x,Jp);
        hold on
        title(name{1},'interpreter','latex','FontSize',16);
        xlabel('Position ($m$)','interpreter','latex','FontSize',14);
        ylabel({'Current Density ($A/m^2$)'},'interpreter','latex','FontSize',14);

        hold off
     end
 end
 
 
