%% code to create the 'TennisCourtMap.jpg' file
%% the TennisBallMap.png file was created with 3 layers in powerpoint
%% 
figure('color','white');
set(gcf,'position',[1 1 5000 2000]);
% draw the court
% https://sportsaspire.com/tennis-court-size
courtXg = [-60 60 60 -60];
courtYg = [-30 -30 30 30];
fill(courtXg,courtYg,[0 0.65 0],'edgecolor','none');
hold on;
courtXw = [-39 39 39 -39];
courtYw = [-18 -18 18 18];
fill(courtXw,courtYw,[1 1 1],'edgecolor','none');
courtXr = [-38.5 38.5 38.5 -38.5];
courtYr = [-13.5 -13.5 -17.5 -17.5];
fill(courtXr,courtYr,[0.65 0.05 0],'edgecolor','none');
fill(courtXr,abs(courtYr),[0.65 0.05 0],'edgecolor','none');
courtXr = [-38.5  -21 -21  -38.5];
courtYr = [-13 -13 13 13];
fill(courtXr,courtYr,[0.65 0.05 0],'edgecolor','none');
fill(abs(courtXr),courtYr,[0.65 0.05 0],'edgecolor','none');
courtXr = [-20.25 -0.5 -0.5 -20.25];
courtYr = [-0.25 -0.25 -13 -13];
fill(courtXr,courtYr,[0.65 0.05 0],'edgecolor','none');
fill(courtXr,abs(courtYr),[0.65 0.05 0],'edgecolor','none');
fill(abs(courtXr),courtYr,[0.65 0.05 0],'edgecolor','none');
fill(abs(courtXr),abs(courtYr),[0.65 0.05 0],'edgecolor','none');
xlim([-60 60]); ylim([-30 30]); set(gca,'xtick',[]); set(gca,'ytick',[])
axis equal
axis tight
