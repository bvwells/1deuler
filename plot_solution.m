clear all;

load 'variables.m';

output_t = variables(1);
numberReports = variables(2);

for report = 1:numberReports 

   figure(report)
   time = (report-1)* (output_t/(numberReports-1));
   
   filename = strcat('Solution', sprintf('%03d',[report]), '.m');

   u = load(filename);
   
   subplot(2,2,1) 
   plot(u(:,1),u(:,2));
   xlabel('x');
   ylabel('\rho');   
   title(strcat('Density (\rho) at time=',num2str(time)));
   
   subplot(2,2,2) 
   plot(u(:,1),u(:,3));
   xlabel('x');
   ylabel('u');   
   title(strcat('Velocity (u) at time=',num2str(time)));

   subplot(2,2,3) 
   plot(u(:,1),u(:,4));
   xlabel('x');
   ylabel('E');
   title(strcat('Pressure (P) at time=',num2str(time)));

   subplot(2,2,4) 
   plot(u(:,1),u(:,5));
   xlabel('x');
   ylabel('E');
   title(strcat('Internal Energy (e) at time=',num2str(time)));
 
 end

figure(numberReports+1)
load 'x_traj.m';
for j=1:100;
   plot(x_traj(:,j),x_traj(:,101));
   hold on;
end
hold off;
axis tight;
xlabel('x');
ylabel('time');