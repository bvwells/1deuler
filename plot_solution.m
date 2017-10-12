clear all;

output_t = 0.2
numberReports = 10

figure(1)
for report = 1:numberReports 
   
   filename = strcat('Solution', sprintf('%03d',[report]), '.m');

   u = load(filename);
    
   plot(u(:,1),u(:,2));

   hold on;
 end

 axis square;
 axis tight;
 xlabel('x');
 ylabel('u');
