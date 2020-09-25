%Allows to plot multiple charge densities
%Running, will bring up a window where can select the files to plot. Takes
%in the input files which are named i.e. 0.77V.txt

clear all

[File,Path]=uigetfile('*.txt','MultiSelect','on');   

N = numel(File);

figure;

h = gobjects(1,N);

for i = 1:N
    
    name= File(1,i);
    str=sprintf('%s', [Path name{1}]);
    format shortG                                              %change formating so doesn't show 0's for e-11 values.
    
    data = importdata(str);
    
    x = data(:,1);  %positions
    p = data(:,3);  %hole densities
    n = data(:,4);  %electron densities
    
    h(2*i-1) = semilogy(x,n, 'LineWidth',1.,'LineStyle', '--','Color','red');
    hold on
    h(2*i) = semilogy(x,p, 'LineWidth',1.,'LineStyle', '--','Color','blue');
    
    
end

matlab.graphics.internal.setPrintPreferences('DefaultPaperPositionMode','manual')
set(gca, 'FontSize', 20) %20 allows to see more labels, on y  axis which makes graph  clearer
xlabel('Position (m)','interpreter','latex','FontSize',26.4);
ylabel({'Hole Density ($m^{-3}$)'},'interpreter','latex','FontSize',26.4);

%legend labels just 2 curve types--> colors for n and p
legend([h(1) h(N)],{'electrons', 'holes'});

