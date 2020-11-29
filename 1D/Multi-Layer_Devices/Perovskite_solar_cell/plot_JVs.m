% allows to plot multiply  JV files on 1 graph
% In the pop-up dialog, select the JV text files that want to plot

clear all

[File,Path]=uigetfile('*.txt','MultiSelect','on');

N = numel(File);  

figure;

 for i = 1:N

    name= File(1,i);
    str=sprintf('%s', [Path name{1}]);   
     format shortG                                              

     data = importdata(str);

     V(:,i) = data(:,1);
     J(:,i) = data(:,2)/10; %the /10 is to convert A/m^2 to mA/cm^2
     
     for j = 1:size(J,1)
         JV(j) = 10*J(j,i)*V(j,i);  %i is defining WHICH JV curve are working with
     end
     [maxPower, I] = min(JV)     % use min--> b/c in solar cell working regime, my J/s are defined as negative, 
                                 % I is the index where the min occurs
                                 
     maxPower = abs(maxPower)
     %max power points                            
     Jmpp = J(I,i)
     Vmpp = V(I,i)
     
     
     %incident power is considered to be 1000 W/m^2 ?     
     PCE = (maxPower/1000)*100 

     plot(V(:,i),J(:,i), 'LineWidth',1.5);
     hold on
 end
    matlab.graphics.internal.setPrintPreferences('DefaultPaperPositionMode','manual') 
     set(gca, 'FontSize', 24)
     xlabel('Voltage (V)','interpreter','latex','FontSize',26.4);
     ylabel({'Current Density ($mA/cm^2$)'},'interpreter','latex','FontSize',26.4);
     axis([0 1.22 -22 2]);
     
 
