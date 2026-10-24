clear all
close all

%[spectrum_name, spectrum_path]=uigetfile('*m12am.erg','Pick one of the spectrum for the station investigated (ex: *am.erg)');    
tic

addpath(pwd);

cd Raw

[file_name, file_path]=uigetfile('*BLANK.csv','Pick blank file in Raw Folder');    
cd(file_path)

root_name=file_name(1:end-22);

% Make a list of all the sample files to be processed
list_cdom=dir([file_path root_name '*Sample.csv']);


% How many samples?
l_cdom=length(list_cdom); 

% List the file names
sname=list_cdom.name  
 


% Import the wavelength vector and pre-allocate a matrix of absrobance
kit=importdata([file_path file_name]); kit=kit.data; kit=flipud(kit);
Wavelength=kit(:,1);

ABSORBANCES=NaN(size(Wavelength,1),length(list_cdom));

clear kit

for i=1:l_cdom
spectrum_file=fullfile(file_path,list_cdom(i).name);    
kit=importdata(spectrum_file); kit=kit.data; kit=flipud(kit);

ABSORBANCES(:,i)=kit(:,2);

end

root_name=list_cdom(i).name(1:end-16);

Wavelength=kit(:,1);

blank_file=fullfile(file_path,file_name);    
blk=importdata(blank_file); blk=blk.data; blk=flipud(blk);
blk=blk(:,2);



%Start the loop for fitting and getting the parameters
for i=1:l_cdom

sampleID=list_cdom(i).name(1:length(list_cdom(i).name)-14);
    
ABSORBANCES(:,i)=ABSORBANCES(:,i)-blk;

%Set range for fit 
range_wave=[500 700];   % set wavelength range in nm low wavelength first, then high wavelength
range=(range_wave(1)-249):(range_wave(2)-249);

%[m,resid,Jac] = nlinfit(Wavelength(101:351),ABSORBANCES(101:351,i),'curfit',[25 .005 -.01]);   % 700 to 450nm
[m,resid,Jac] = nlinfit(Wavelength(range),ABSORBANCES(range,i),'curfit_CDOM',[25 .005 0.01]); % 700 to 400nm
% curfit is another .m file (curfit.m) which includes the formula for the
% model used (here an exponential with three parameters,... y=a*exp(-b*x)-c where a, b and c are parameters and x and y are variables)


%  Calculate yhat - fitted abs
[yhat] = curfit_CDOM(m,Wavelength);
fitted_abs(:,i)=yhat;
% Subtract paramter 3 which is the offset value

ABS(:,i)=ABSORBANCES(:,i)-m(3);


%pathlength in m
pth=0.01;

% Put all offset values in one vector
nnn(i)=m(3);

% Convert absobances to absorption coefficients
ag(:,i)=ABS(:,i).*log(10)./pth;

end

% End of loop




% Set colors for the plotting
colors=[1 0 0 ; 1 0.55 0 ; 1 1 0 ; 0 0 1 ; 0 1 1 ; 0 1 0 ; 1 0 1];

if l_cdom >= 8    
    for i=8:l_cdom
        colors(i,:)= [rand([1 1]) rand([1 1]) rand([1 1])]  
    end
else
 colors=colors;   
end    

% Get the names of samples right
    for i=1:l_cdom
        list_cdom(i).name=list_cdom(i).name(1:length(list_cdom(i).name)-11);
    end



%% Start ploting

figure('Position',[20,50,1600,1200])
subplot(2,3,1)
hold
for i=1:l_cdom
    plot(Wavelength(51:351),(ABSORBANCES(51:351,i)),'-','Color',colors(i,:),'LineWidth',2)   
end
for i=1:l_cdom
    plot(Wavelength(51:351),(fitted_abs(51:351,i)),'-.','Color',colors(i,:)*0.8,'LineWidth',2)   
end
xlabel('Wavelength (nm)','FontSize',14)
ylabel('Absorbance (unitless)','FontSize',14)
axis tight
title('Absorbances before baseline adjustment','FontSize',14)
h=legend(list_cdom(1:length(list_cdom)).name);
set(h,'Interpreter','none')
grid on


subplot(2,3,2)
hold
for i=1:l_cdom
    plot(Wavelength(251:451),(ABSORBANCES(251:451,i)),'-','Color',colors(i,:),'LineWidth',2)   
end
for i=1:l_cdom
    plot(Wavelength(251:451),(fitted_abs(251:451,i)),'-.','Color',colors(i,:)*0.8,'LineWidth',2)   
end
xlabel('Wavelength (nm)','FontSize',14)
ylabel('Absorbance (unitless)','FontSize',14)
axis tight
title('Absorbances before baseline adjustment - Long-\lambda range','FontSize',14)
h=legend(list_cdom(1:length(list_cdom)).name);
set(h,'Interpreter','none')
grid on

subplot(2,3,3)
hold
for i=1:l_cdom
    plot(Wavelength(1:451),(ag(1:451,i)),'-','Color',colors(i,:),'LineWidth',2) 
end
xlabel('Wavelength (nm)','FontSize',14)
ylabel('Absorption coefficient of CDOM (m^-^1)','FontSize',14)
axis tight
title('Absorption Coefficient - Full spectral range','FontSize',14)
h=legend(list_cdom(1:length(list_cdom)).name);
set(h,'Interpreter','none')
grid on


subplot(2,3,4)
hold
for i=1:l_cdom
    plot(Wavelength(1:151),(ag(1:151,i)),'-','Color',colors(i,:),'LineWidth',2) 
end
xlabel('Wavelength (nm)','FontSize',14)
ylabel('Absorption coefficient of CDOM (m^-^1)','FontSize',14)
axis tight
title('Absorption Coefficient - UV range','FontSize',14)
h=legend(list_cdom(1:length(list_cdom)).name);
set(h,'Interpreter','none')
grid on


subplot(2,3,5)
hold
for i=1:l_cdom
    plot(Wavelength(51:251),(ag(51:251,i)),'-','Color',colors(i,:),'LineWidth',2) 
end
xlabel('Wavelength (nm)','FontSize',14)
ylabel('Absorption coefficient of CDOM (m^-^1)','FontSize',14)
axis tight
title('Absorption Coefficient - 350-550-nm range','FontSize',14)
h=legend(list_cdom(1:length(list_cdom)).name);
set(h,'Interpreter','none')
grid on


subplot(2,3,6)
hold
for i=1:l_cdom
    plot(Wavelength(251:551),(ag(251:551,i)),'-','Color',colors(i,:),'LineWidth',2) 
end
xlabel('Wavelength (nm)','FontSize',14)
ylabel('Absorption coefficient of CDOM (m^-^1)','FontSize',14)
axis tight
title('Absorption Coefficient - Long-\lambda range','FontSize',14)
h=legend(list_cdom(1:length(list_cdom)).name);
set(h,'Interpreter','none')
grid on


cops_ind=[305 320 340 380 395 412 443 465 490 510 532 560 589 625 665 683 694 710 780]-249;

% Saving the data
ag_wave=[Wavelength ag];
ag_wave_cops=[Wavelength(cops_ind) ag(cops_ind,:)];



labels_columns=cell(1,l_cdom);
labels_columns(1,1)=cellstr('Wavelength');
for i=2:l_cdom+1
        labels_columns(1,i)=[cellstr(list_cdom(i-1).name)];   
end





% Save the file
Dataset=[[cellstr(labels_columns)]' num2cell([ag_wave(1:551,:)'])];
Dataset_cops=[[cellstr(labels_columns)]' num2cell([ag_wave_cops'])];


for i=1:l_cdom

labels_s=cell(1,4);
labels_s(1,1)=cellstr('S275295');
labels_s(1,2)=cellstr('S350400');
labels_s(1,3)=cellstr('Sr');
labels_s(1,4)=cellstr('S300400');

range=[300:400];
range1=[275:295];
range2=[350:400];

S=robustfit(Wavelength(range-249),log(ag(range-249,i)));
S300400=S(2);

S1=robustfit(Wavelength(range1-249),log(ag(range1-249,i)));
S275295=S1(2);

S2=robustfit(Wavelength(range2-249),log(ag(range2-249,i)));
S350400=S2(2);

Sr=S275295./S350400;

Slopes(i,:)=[-S275295 -S350400 Sr -S300400]; 

end

Slopes2=[labels_s ;num2cell([Slopes])]; 


labels_columns=cell(1,l_cdom);
labels_columns(1,1)=cellstr('Slopes');
for i=2:l_cdom+1
        labels_columns(1,i)=[cellstr(list_cdom(i-1).name)];   
end



Dataset_slopes=[[cellstr(labels_columns)]' Slopes2];




% Move in Processed Folder
cd ../Processed/

% Save the output in a new textefile
save_name=['CDOM_fitted_offset_' root_name(1:end-1) '.csv'];
dlmcell(save_name,Dataset,',')

save_name=['CDOM_fitted_offset_' root_name(1:end-1) '_for_C-OPS.csv'];
dlmcell(save_name,Dataset_cops,',')

save_name=['CDOM_fitted_offset_' root_name(1:end-1) '_SLOPES.csv'];
dlmcell(save_name,Dataset_slopes,',')

cd ..