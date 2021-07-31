%% Test against Imre's time dependent exact solution in 1D, 
%% Original CNe, hop and leap-hop combinations 
clear;close all;
clear; fprintf('\n')
tic
Mh=18; %How many different timestep we want 
D=1; x0=-5; xf=5; ti=0.5; tf=1; TIME=tf-ti; %Diffusion coefficient, lenght of the space and time interval, 
N=1001; Dx=(xf-x0)/(N-1);   % number of cells, Space step size
MaxD=zeros(Mh,7);Axhstep=zeros(Mh,1); rfun=zeros(Mh,1); %%Errors and axis for the timestep-sizes
T=5; %% initial number of time levels
xaxis=zeros(N,1); for i=1:N xaxis(i)=x0+(i-1)*Dx; end  %x axis
M1=zeros(N,N); ee1=zeros(N,1);ee2=zeros(N,1);
u0=zeros(N,1); uex=zeros(N,1); %%initial values and exact solution
E=D/(Dx)^2; %%matrix elements
k=zeros(N,1); %%sources, now zero
cpar=zeros(1,N);  %%cell parity vector
for i=1:N 
 u0(i)=xaxis(i)/ti^(2)*exp(-xaxis(i)^2/2/D/ti^2)*hypergeom(1/2,3/2,xaxis(i)^2/2/D/ti^2); 
 uex(i)=xaxis(i)/tf^(2)*exp(-xaxis(i)^2/2/D/tf^2)*hypergeom(1/2,3/2,xaxis(i)^2/2/D/tf^2);   %%Initial and exact final results
 %u0(i)=xaxis(i)/ti^(2)*exp(-xaxis(i)^2/2/D/ti^2)*kummerU(1/2,3/2,xaxis(i)^2/2/D/ti^2); 
 %uex(i)=xaxis(i)/tf^(2)*exp(-xaxis(i)^2/2/D/tf^2)*kummerU(1/2,3/2,xaxis(i)^2/2/D/tf^2);
end
plot(xaxis,u0,xaxis,uex)
M1(1,1)=-E;M1(1,2)=E; M1(N,N)=-E;M1(N,N-1)=E;
for i=2:N-1     M1(i,i-1)=E; M1(i,i)=-2*E; M1(i,i+1)=E; end
Neb1=zeros(N,1);Neb2=zeros(N,1); 
for i=1:N
    if (i>1) Neb1(i)=M1(i,i-1); end  %bal szomszéd
    if (i<N) Neb2(i)= M1(i,i+1); end %jobb szomszéd
   end
%%%%%%%%%%%%
EG = eig(M1);Emax=0; E0=-1; Emin=-2; %%Eigenvalues
for i=1:N if EG(i)<Emax Emax=EG(i); end
    if EG(i)>E0 
       E0=EG(i); end 
end
for i=1:N if EG(i)>Emin
      if EG(i)<E0 
           Emin=EG(i); end 
    end
end
hMAX= -2/Emax  %% max timestep for Expl Euler
StiffRatio=Emax/Emin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ih=1:Mh  %% Big loop for the timestep, ih is the Serial Number of run, 1 for the largest timestep; 
ih
h=TIME/(T-1);
%    T = ceil(TIME/h-0.1)+1;
tax=zeros(T,1); rt=zeros(T,1);leftbord=zeros(T,1);rightbord=zeros(T,1); %%physical time
tic
for t=1:T
      tax(t)=ti+(t-1)*h;
      leftbord(t)=x0/tax(t)^(2)*exp(-x0^2/2/D/tax(t)^2)*hypergeom(1/2,3/2,x0^2/2/D/tax(t)^2);
      rightbord(t)=xf/tax(t)^(2)*exp(-xf^2/2/D/tax(t)^2)*hypergeom(1/2,3/2,xf^2/2/D/tax(t)^2);
      rt(t)=h*E*(tax(t)+h/2);
      eet(t)=exp(-2*rt(t)); 
      ee2t(t)=exp(-rt(t));
end
toc
if(h>0.00001 && h<0.0001) rfun(ih)=2000000000*h*h*h; end
UC=zeros(N,1); UTEMP1=zeros(N,1); UTEMP2=zeros(N,1);  %%Final U and temporary results
UL1=zeros(N,1);UL2=zeros(N,1);UL3=zeros(N,1);UL4=zeros(N,1);UL5=zeros(N,1);
a0=zeros(N,1);a1=zeros(N,1);s1=zeros(N,1);
Axhstep(ih)=h; %%Timestep axis
r=h*E; ee2=ee2t(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UC(:)=u0(:); %%%% 1 Stage CNe
tic
for t=1:T-1
    for i=2:N-1  %%%Neighbours are constans
        a0(i)=(UC(i-1)+UC(i+1))/2;
        UTEMP1(i)=UC(i)*eet(t)+ a0(i)*(1-eet(t));
    end
    UTEMP1(1)=leftbord(t+1); 
    UTEMP1(N)=rightbord(t+1);
 UC(:)=UTEMP1(:);
end
toc
%%Constant errors
Kul=zeros(N,1); MaxKul=0; 
for i=1:N
     Kul(i)=UC(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,1)=MaxKul;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% L1 L1
UL1=u0(:); %%%% A 20 method Leap Hop, full CNe
tic
 for i=2:2:N-1  %%% even
      UL1(i)=UL1(i)*ee2+ ((UL1(i-1)+UL1(i+1))/2)*(1-ee2);
  end
for t=1:2:T-2
    for i=3:2:N-2  %%% odd
    UL1(i)=UL1(i)*eet(t+1)+ ((UL1(i-1)+UL1(i+1))/2)*(1-eet(t+1)); 
    end
    UL1(1)=leftbord(t+1);   UL1(N)=rightbord(t+1);
    for i=2:2:N-1  %%% even
         UL1(i)=UL1(i)*eet(t+1)+ ((UL1(i-1)+UL1(i+1))/2)*(1-eet(t+1));    
    end
    for i=3:2:N-2  %%% odd
    UL1(i)=UL1(i)*eet(t+2)+ ((UL1(i-1)+UL1(i+1))/2)*(1-eet(t+2));    
    end
    UL1(1)=leftbord(t+2);     UL1(N)=rightbord(t+2);    
    the=1/2;
    for i=2:2:N-1  %%% even
             if (tax(t+1)<tf-h)
    UL1(i)=UL1(i)*eet(t+2)+ ((UL1(i-1)+UL1(i+1))/2)*(1-eet(t+2));
            else
    UL1(i)=UL1(i)*ee2t(t+2)+ ((UL1(i-1)+UL1(i+1))/2)*(1-ee2t(t+2));
            end
    end
 end
toc
%% C5 errors
Kul=zeros(N,1); MaxKul=0; 
for i=1:N
     Kul(i)=UL1(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,2)=MaxKul;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% L2, L2
UL2=u0(:); %%%% L2 method Leap Hop, (0, 1/2, 1/2, 1/2, 1/2) 
tic
the=0;
 for i=2:2:N-1  %%% UPFD, even
      UL2(i)=((1-rt(1)*the)*UL2(i)+rt(1)/2*(UL2(i-1)+UL2(i+1)))/(1+rt(1)*(1-the));
 end
for t=1:2:T-2
    the=1/2;
    for i=3:2:N-2  %%% theta=1/2, odd
    UL2(i)=((1-2*rt(t+1)*the)*UL2(i)+rt(t+1)*(UL2(i-1)+UL2(i+1)))/(1+2*rt(t+1)*(1-the));    
    end
    UL2(1)=leftbord(t+1); 
    UL2(N)=rightbord(t+1);
    the=1/2;
    for i=2:2:N-1  %%% theta=1/2, even
    UL2(i)=((1-2*rt(t+1)*the)*UL2(i)+rt(t+1)*(UL2(i-1)+UL2(i+1)))/(1+2*rt(t+1)*(1-the));
    end
    the=1/2;
    for i=3:2:N-2  %%% 
    UL2(i)=((1-2*rt(t+2)*the)*UL2(i)+rt(t+2)*(UL2(i-1)+UL2(i+1)))/(1+2*rt(t+2)*(1-the)); 
    end
    UL2(1)=leftbord(t+2); 
    UL2(N)=rightbord(t+2);    
    the=1/2;
    for i=2:2:N-1  %%% 
             if (tax(t+1)<tf-h)
    UL2(i)=((1-2*rt(t+2)*the)*UL2(i)+rt(t+2)*(UL2(i-1)+UL2(i+1)))/(1+2*rt(t+2)*(1-the));
            else
    UL2(i)=((1-rt(t+2)*the)*UL2(i)+rt(t+2)/2*(UL2(i-1)+UL2(i+1)))/(1+rt(t+2)*(1-the));
            end
    end
end
toc
%% C5 errors
Kul=zeros(N,1); MaxKul=0; 
for i=1:N
     Kul(i)=UL2(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,3)=MaxKul;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% L3 L3
UL3=u0(:); %%%% L3 method Leap Hop (1/5, 1/2, 1/2, 1/2, 1/2) 
tic
the=1/5;
 for i=2:2:N-1  %%% UPFD
       UL3(i)=((1-rt(1)*the)*UL3(i)+rt(1)/2*(UL3(i-1)+UL3(i+1)))/(1+rt(1)*(1-the));
 end
for t=1:2:T-2
    the=1/2;
    for i=3:2:N-2  %%%
    UL3(i)=((1-2*rt(t+1)*the)*UL3(i)+rt(t+1)*(UL3(i-1)+UL3(i+1)))/(1+2*rt(t+1)*(1-the));    
    end
    UL3(1)=leftbord(t+1); 
    UL3(N)=rightbord(t+1);
    the=1/2;
    for i=2:2:N-1  %%%
    UL3(i)=((1-2*rt(t+1)*the)*UL3(i)+rt(t+1)*(UL3(i-1)+UL3(i+1)))/(1+2*rt(t+1)*(1-the));
    end
    the=1/2;
    for i=3:2:N-2  %%% 
    UL3(i)=((1-2*rt(t+2)*the)*UL3(i)+rt(t+2)*(UL3(i-1)+UL3(i+1)))/(1+2*rt(t+2)*(1-the)); 
    end
    UL3(1)=leftbord(t+2); 
    UL3(N)=rightbord(t+2);    
    the=1/2;
    for i=2:2:N-1  %%% 
             if (tax(t+1)<tf-h)
    UL3(i)=((1-2*rt(t+2)*the)*UL3(i)+rt(t+2)*(UL3(i-1)+UL3(i+1)))/(1+2*rt(t+2)*(1-the));
            else
    UL3(i)=((1-rt(t+2)*the)*UL3(i)+rt(t+2)/2*(UL3(i-1)+UL3(i+1)))/(1+rt(t+2)*(1-the));
            end
    end
end
toc
%% C5 errors
Kul=zeros(N,1); MaxKul=0; 
for i=1:N
     Kul(i)=UL3(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,4)=MaxKul;
%%%%%%%%%%%%%%%%%%%%%%%%% L4 L4
UL4=u0(:); %%%% L4 method Leap Hop (1/4, 1/2, CNe, 1/2, 1/2) 
tic
the=1/4;
 for i=2:2:N-1  %%% UPFD
     UL4(i)=((1-rt(1)*the)*UL4(i)+rt(1)/2*(UL4(i-1)+UL4(i+1)))/(1+rt(1)*(1-the));
 end
for t=1:2:T-2
    the=1/2;
    for i=3:2:N-2  %%% CrN
    UL4(i)=((1-2*rt(t+1)*the)*UL4(i)+rt(t+1)*(UL4(i-1)+UL4(i+1)))/(1+2*rt(t+1)*(1-the));    
    end
    UL4(1)=leftbord(t+1); 
    UL4(N)=rightbord(t+1);
    for i=2:2:N-1  %%% theta=1/3
         UL4(i)=UL4(i)*eet(t+1)+ ((UL4(i-1)+UL4(i+1))/2)*(1-eet(t+1));    
    end
    the=1/2;
    for i=3:2:N-2  %%% CNe
    UL4(i)=((1-2*rt(t+2)*the)*UL4(i)+rt(t+2)*(UL4(i-1)+UL4(i+1)))/(1+2*rt(t+2)*(1-the)); 
    end
    UL4(1)=leftbord(t+2); 
    UL4(N)=rightbord(t+2);    
    the=1/2;
    for i=2:2:N-1  %%% CrN
             if (tax(t+1)<tf-h)
    UL4(i)=((1-2*rt(t+2)*the)*UL4(i)+rt(t+2)*(UL4(i-1)+UL4(i+1)))/(1+2*rt(t+2)*(1-the));
            else
    UL4(i)=((1-rt(t+2)*the)*UL4(i)+rt(t+2)/2*(UL4(i-1)+UL4(i+1)))/(1+rt(t+2)*(1-the));
            end
    end
end
toc
%% C5 errors
Kul=zeros(N,1); MaxKul=0; 
for i=1:N
     Kul(i)=UL4(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,5)=MaxKul;
%%%%%%%%%%%%%%%%%%%%%%%%% L5 L5
UL5=u0(:); %%%% L4 method Leap Hop (1/5, 1/2, CNe, 1/2, 1/2) 
tic
the=1/5;
 for i=2:2:N-1  %%% 
    UL5(i)=((1-rt(1)*the)*UL5(i)+rt(1)/2*(UL5(i-1)+UL5(i+1)))/(1+rt(1)*(1-the));
 end
for t=1:2:T-2
    the=1/2;
    for i=3:2:N-2  %%% 
    UL5(i)=((1-2*rt(t+1)*the)*UL5(i)+rt(t+1)*(UL5(i-1)+UL5(i+1)))/(1+2*rt(t+1)*(1-the));    
    end
    UL5(1)=leftbord(t+1); 
    UL5(N)=rightbord(t+1);
    for i=2:2:N-1  %%% CNe
         UL5(i)=UL5(i)*eet(t+1)+ ((UL5(i-1)+UL5(i+1))/2)*(1-eet(t+1));    
      end
    the=1/2;
    for i=3:2:N-2  %%% 
    UL5(i)=((1-2*rt(t+2)*the)*UL5(i)+rt(t+2)*(UL5(i-1)+UL5(i+1)))/(1+2*rt(t+2)*(1-the)); 
    end
    UL5(1)=leftbord(t+2); 
    UL5(N)=rightbord(t+2);    
    the=1/2;
    for i=2:2:N-1  %%% 
             if (tax(t+1)<tf-h)
    UL5(i)=((1-2*rt(t+2)*the)*UL5(i)+rt(t+2)*(UL5(i-1)+UL5(i+1)))/(1+2*rt(t+2)*(1-the));
            else
    UL5(i)=((1-rt(t+2)*the)*UL5(i)+rt(t+2)/2*(UL5(i-1)+UL5(i+1)))/(1+rt(t+2)*(1-the));
            end
    end
end
toc
%% C5 errors
Kul=zeros(N,1); MaxKul=0; 
for i=1:N
     Kul(i)=UL5(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,6)=MaxKul;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% A random combination to test
% UC5=u0(:); %%%% C5 method Leap Hop
% tic
% the=0.1;
%  for i=2:2:N-1  %%% UPFD
%       UC5(i)=UC5(i)*ee2+ ((UC5(i-1)+UC5(i+1))/2)*(1-ee2);
%       %UC5(i)=((1-rt(t)*the)*UC5(i)+rt(t)/2*(UC5(i-1)+UC5(i+1)))/(1+rt(t)*(1-the));
%  end
% for t=1:2:T-2
%     the=0.1;
%     for i=3:2:N-2  %%% CrN
%     UC5(i)=((1-2*rt(t+1)*the)*UC5(i)+rt(t+1)*(UC5(i-1)+UC5(i+1)))/(1+2*rt(t+1)*(1-the));    
%     %UC5(i)=(r/(1+r))*(UC5(i-1)+UC5(i+1))+((1-r)/(1+r))*UC5(i);
%     end
%     UC5(1)=leftbord(t+1); 
%     UC5(N)=rightbord(t+1);
%     the=0.4;
%     for i=2:2:N-1  %%% theta=1/3
%          %UC5(i)=UC5(i)*eet(t)+ ((UC5(i-1)+UC5(i+1))/2)*(1-eet(t));    
%         UC5(i)=((1-2*rt(t+1)*the)*UC5(i)+rt(t+1)*(UC5(i-1)+UC5(i+1)))/(1+2*rt(t+1)*(1-the));
%     end
%     the=1/2;
%     for i=3:2:N-2  %%% CNe
%     %UC5(i)=UC5(i)*eet(t+2)+ ((UC5(i-1)+UC5(i+1))/2)*(1-eet(t+2));    
%     UC5(i)=((1-2*rt(t)*the)*UC5(i)+rt(t)*(UC5(i-1)+UC5(i+1)))/(1+2*rt(t)*(1-the)); 
%     end
%     UC5(1)=leftbord(t+2); 
%     UC5(N)=rightbord(t+2);    
%     the=0.3;
%     for i=2:2:N-1  %%% CrN
%              if (tax(t+1)<tf-h)
%     UC5(i)=((1-2*rt(t+2)*the)*UC5(i)+rt(t+2)*(UC5(i-1)+UC5(i+1)))/(1+2*rt(t+2)*(1-the));
%             else
%     UC5(i)=((1-rt(t+2)*the)*UC5(i)+rt(t+2)/2*(UC5(i-1)+UC5(i+1)))/(1+rt(t+2)*(1-the));
%             end
%     end
% end
% toc
% %% C5 errors
% Kul=zeros(N,1); MaxKul=0; 
% for i=1:N
%      Kul(i)=UC5(i)-uex(i);
%     if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
% end
% MaxD(ih,7)=MaxKul;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%88888888888888888888888


%%%%%%%%%%%%%%%%%%%%%%%%%%
%h=h/fact;
T=2*(T-1)+1;
end

    green         = [0 1 0.4]; 
    dark_grey     = [0.5,0.5,0.5]; 
    brown         = [0.2 0 0];
    orange        = [1 0.5 0];
    dark_red      = [0.75, 0.0780, 0.1840]; 
    dark_yellow   = [0.9290, 0.6940, 0.1250];
    yellow   = [0.990, 0.8940, 0.250];  
    dark_green    = [0, 0.45, 0];
    dark_purple   = [0.55, 0.0840, 0.5560];
    light_blue     = [0.3, 0.65, 1];
    dark_dark_grey = [0.2, 0.2, 0.25]; 

Magenta	 =    [1, 0, 1];
MilkChocolate	= [165/255, 42/255, 42/255];
Black		= [0, 0, 0];
Red	= [1, 0, 0];
GoldenBrown =	[0.9, 0.75, 0.1];
Blue	=  [0, 0, 1]; 
Green = [0/255, 128/255, 0/255]; 
%Green		#008000		0.33, 1.00, 0.25];

plot(xaxis,u0, '--',  'Color', dark_grey, 'LineWidth', 1);
	hold on;
plot(xaxis,uex, '-',  'Color', dark_red, 'LineWidth', 2);
	hold on;
plot(xaxis,UL2, ':',  'Color', green, 'LineWidth', 6);
	hold on;
hold off;
XLabel = xlabel('x');
set(XLabel, 'FontSize', 14);
set(XLabel, 'FontWeight', 'bold');

YLabel = ylabel('u');
set(YLabel, 'FontSize', 14); 
set(YLabel, 'FontWeight', 'bold');

Legend=legend({'u0','u^{ref}', 'L2'},'Location','southeast', 'FontSize', 13);
	%% figure 1:
  figure('Name', 'Errors as a fuction of time step h');   

	plot(Axhstep(:),MaxD(:,1), ":o", 'Color',Magenta, 'LineWidth', 1.8, 'MarkerSize',5); %% cn 
	hold on;
    plot(Axhstep(:),MaxD(:,2), ':',  'Color', dark_grey, 'LineWidth', 5); %% cncn p=1/3
	hold on;
    plot(Axhstep(:),MaxD(:,3), '-x',  'Color', light_blue, 'LineWidth', 2.7); %% cncn p=2/3
	hold on;
    plot(Axhstep(:),MaxD(:,4), ':O',  'Color', dark_yellow, 'LineWidth', 3); %% cncn p=1/2
	hold on;
	plot(Axhstep(:),MaxD(:,5), '-bd',  'Color', MilkChocolate, 'LineWidth', 1.1,'MarkerSize',6); %% cncn p=3/4
	hold on;
    plot(Axhstep(:),MaxD(:,6), '-.v',  'Color', Black, 'LineWidth', 1.0, 'MarkerSize',4); %% cncn p=3/4
	hold on;

	hold off;
		%% Here pl=6-> MaxInndex=pl+2=6+2=8 . 0_Duffort-Frankel 


XLabel = xlabel('Effective time step size h_{EFF}');
set(XLabel, 'FontSize', 16);
set(XLabel, 'FontWeight', 'bold');

YLabel = ylabel('Errors');
set(YLabel, 'FontSize', 16); 
set(YLabel, 'FontWeight', 'bold');


% Get handle to current axes.
Ax = gca;
Ax.XAxis.FontSize = 18;
Ax.YAxis.FontSize = 18;
Ax.FontWeight = 'bold';
Ax.TickLength = [0.018 0.035];
Ax.LineWidth = 1.2;

Legend=legend({'CNe','L1', 'L2','L3','L4','L5'},'Location','southeast', 'FontSize', 14);
%Legend=legend({'CNe','A1 (C;C;C;C;C) ','A2 (1/4;1/2;C;1/2;3/4)','A3 (1/4;1/2;1/2;1/2;3/4)','A4 (0;1/2;1/2;1/2;1)','A5 (0;1/2;1/2;C;1)'},'Location','southeast');
set(Legend, 'FontSize', 14); 
set(Legend, 'FontWeight', 'bold');
%set(Legend, 'NumColumns', 2);
%title(Legend, 'My Legend Title');
set(gca,'xScale','log');
set(gca,'yScale','log');
grid on;

