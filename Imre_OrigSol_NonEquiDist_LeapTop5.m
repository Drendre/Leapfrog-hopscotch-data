%% Solves the diffusion equation, 1D, space-dependent grid, LEAPHOP + CNe + Heun 
%% ONLY After running the Generating non-equidistant grid code 
Mh=18;
MaxD=zeros(Mh,8);Axhstep=zeros(Mh,1); rfun=zeros(Mh,1);%%Errors and axis for the timestep-sizes
T=5;
cpar=zeros(1,N);  %%cell parity vector
Neb1=zeros(N,1);Neb2=zeros(N,1); 
ee=zeros(N,1);ee2=zeros(N,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ih=1:Mh  %% Big loop for the timestep, ih is the Serial Number of run, 1 for the largest timestep; 
ih
h=TIME/(T-1); Axhstep(ih)=h; %%Timestep axis
for i=1:N
    if (i>1) Neb1(i)=h*M1(i,i-1); end  %left neighbor
    if (i<N) Neb2(i)=h*M1(i,i+1); end %right neighbor
end
tax=zeros(T,1); %%physical time axis
rN=zeros(N,1); %% r depends on space
for t=1:T
      tax(t)=ti+(t-1)*h; 
end
 for i=1:N
     rN(i)=-h*M1(i,i);
     ee(i)=exp(-rN(i)); 
     ee2(i)=exp(-rN(i)/2);
  end
%if(h>0.00001 && h<0.0001) rfun(ih)=2000000000*h*h; end
UC=zeros(N,1);UTr=zeros(N,1);UTEMP=zeros(N,1);UTEMP2=zeros(N,1);
UL1=zeros(N,1);UL2=zeros(N,1);UL3=zeros(N,1);UL4=zeros(N,1);UL5=zeros(N,1);
UC(:)=u0(:);UL1(:)=u0(:);UL2(:)=u0(:);UL3(:)=u0(:);UL4(:)=u0(:);UL5(:)=u0(:); %% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UTr(:)=u0(:); %%%% 2 stage explicit trapezoidal, Heun
tic
for t=1:T-1
     A=zeros(N,1);
    for i=1:N  %%% Exp Euler
       if(i>1) A(i)=A(i)+Neb1(i)*UTr(i-1); end
       if(i<N) A(i)=A(i)+Neb2(i)*UTr(i+1); end
        UTEMP(i)=(1-rN(i))*UTr(i)+A(i);
    end
     A=zeros(N,1);
    for i=1:N  %%% Trapezoidal
      if(i>1) A(i)=A(i)+Neb1(i)*(UTr(i-1)+UTEMP(i-1))/2; end
      if(i<N) A(i)=A(i)+Neb2(i)*(UTr(i+1)+UTEMP(i+1))/2; end    
   UTEMP2(i)=UTr(i)-rN(i)*(UTr(i)+UTEMP(i))/2 + A(i);
    end
UTr=UTEMP2;
end
toc
%%Heun errors
Kul=zeros(N,1); MaxKul=0; 
for i=1:N
     Kul(i)=UTr(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,7)=MaxKul; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
 for t=1:T-1
     A=zeros(N,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CNe ORIG
    for i=1:N  %%%Neighbours are constans
       if(i>1) A(i)=A(i)+Neb1(i)*UC(i-1); end
       if(i<N) A(i)=A(i)+Neb2(i)*UC(i+1); end
           UTEMP(i)=UC(i)*ee(i)+ A(i)/rN(i)*(1-ee(i));
    end
 UC(:)=UTEMP(:);
 end
%% CNe errors
Kul=zeros(N,1); MaxKul=0; 
for i=1:N
     Kul(i)=UC(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,1)=MaxKul;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% L1 L1
  for i=1:2:N  %%% odd
    A=zeros(N,1);
      if(i>1) A(i)=A(i)+Neb1(i)*UL1(i-1); end
      if(i<N) A(i)=A(i)+Neb2(i)*UL1(i+1); end
      UL1(i)=UL1(i)*ee2(i)+ A(i)/rN(i)*(1-ee2(i)); % CNe
    A=zeros(N,1); the=0;
      if(i>1) A(i)=A(i)+Neb1(i)*UL2(i-1); end
      if(i<N) A(i)=A(i)+Neb2(i)*UL2(i+1); end
      UL2(i)=((1-rN(i)/2*the)*UL2(i)+ A(i)/2)/(1+rN(i)/2*(1-the)); % Theta 0
    A=zeros(N,1); the=1/5;
      if(i>1) A(i)=A(i)+Neb1(i)*UL3(i-1); end
      if(i<N) A(i)=A(i)+Neb2(i)*UL3(i+1); end
      UL3(i)=((1-rN(i)/2*the)*UL3(i)+ A(i)/2)/(1+rN(i)/2*(1-the)); % Theta
    A=zeros(N,1); the=1/4;
      if(i>1) A(i)=A(i)+Neb1(i)*UL4(i-1); end
      if(i<N) A(i)=A(i)+Neb2(i)*UL4(i+1); end
      UL4(i)=((1-rN(i)/2*the)*UL4(i)+ A(i)/2)/(1+rN(i)/2*(1-the)); % Theta
    A=zeros(N,1); the=1/5;
      if(i>1) A(i)=A(i)+Neb1(i)*UL5(i-1); end
      if(i<N) A(i)=A(i)+Neb2(i)*UL5(i+1); end
      UL5(i)=((1-rN(i)/2*the)*UL5(i)+ A(i)/2)/(1+rN(i)/2*(1-the)); % Theta
 
 end
  %%%%%%%%%%%%%%%% BIG LOOP FOR TIME FOR THE LEAPFROG 
for t=1:2:T-2
  
   for i=2:2:N  %%%  even
    A=zeros(N,1);
      if(i>1) A(i)=A(i)+Neb1(i)*UL1(i-1); end
      if(i<N) A(i)=A(i)+Neb2(i)*UL1(i+1); end
      UL1(i)=UL1(i)*ee(i)+ A(i)/rN(i)*(1-ee(i)); % CNe
    A=zeros(N,1); the=1/2;
      if(i>1) A(i)=A(i)+Neb1(i)*UL2(i-1); end
      if(i<N) A(i)=A(i)+Neb2(i)*UL2(i+1); end
     UL2(i)=((1-rN(i)*the)*UL2(i)+ A(i))/(1+rN(i)*(1-the)); % Theta 1/2
     A=zeros(N,1); the=1/2;
      if(i>1) A(i)=A(i)+Neb1(i)*UL3(i-1); end
      if(i<N) A(i)=A(i)+Neb2(i)*UL3(i+1); end
     UL3(i)=((1-rN(i)*the)*UL3(i)+ A(i))/(1+rN(i)*(1-the)); % Theta
     A=zeros(N,1); the=1/2;
      if(i>1) A(i)=A(i)+Neb1(i)*UL4(i-1); end
      if(i<N) A(i)=A(i)+Neb2(i)*UL4(i+1); end
     UL4(i)=((1-rN(i)*the)*UL4(i)+ A(i))/(1+rN(i)*(1-the)); % Theta
      A=zeros(N,1); the=1/2;
      if(i>1) A(i)=A(i)+Neb1(i)*UL5(i-1); end
      if(i<N) A(i)=A(i)+Neb2(i)*UL5(i+1); end
     UL5(i)=((1-rN(i)*the)*UL5(i)+ A(i))/(1+rN(i)*(1-the)); % Theta
   end
    for i=1:2:N  %%% CNe odd
     A=zeros(N,1);
      if(i>1) A(i)=A(i)+Neb1(i)*UL1(i-1); end
      if(i<N) A(i)=A(i)+Neb2(i)*UL1(i+1); end
      UL1(i)=UL1(i)*ee(i)+ A(i)/rN(i)*(1-ee(i));
       A=zeros(N,1); the=1/2;
      if(i>1) A(i)=A(i)+Neb1(i)*UL2(i-1); end
      if(i<N) A(i)=A(i)+Neb2(i)*UL2(i+1); end
     UL2(i)=((1-rN(i)*the)*UL2(i)+ A(i))/(1+rN(i)*(1-the)); % Theta 1/2
     A=zeros(N,1); the=1/2;
      if(i>1) A(i)=A(i)+Neb1(i)*UL3(i-1); end
      if(i<N) A(i)=A(i)+Neb2(i)*UL3(i+1); end
     UL3(i)=((1-rN(i)*the)*UL3(i)+ A(i))/(1+rN(i)*(1-the)); % Theta
     A=zeros(N,1); 
      if(i>1) A(i)=A(i)+Neb1(i)*UL4(i-1); end
      if(i<N) A(i)=A(i)+Neb2(i)*UL4(i+1); end
     UL4(i)=UL4(i)*ee(i)+ A(i)/rN(i)*(1-ee(i)); % CNe
      A=zeros(N,1); the=1/2;
      if(i>1) A(i)=A(i)+Neb1(i)*UL5(i-1); end
      if(i<N) A(i)=A(i)+Neb2(i)*UL5(i+1); end
     UL5(i)=UL5(i)*ee(i)+ A(i)/rN(i)*(1-ee(i)); % CNe
    end
  for i=2:2:N  %%% CNe even
      A=zeros(N,1);
      if(i>1) A(i)=A(i)+Neb1(i)*UL1(i-1); end
      if(i<N) A(i)=A(i)+Neb2(i)*UL1(i+1); end
      UL1(i)=UL1(i)*ee(i)+ A(i)/rN(i)*(1-ee(i)); % CNe
    A=zeros(N,1); the=1/2;
      if(i>1) A(i)=A(i)+Neb1(i)*UL2(i-1); end
      if(i<N) A(i)=A(i)+Neb2(i)*UL2(i+1); end
     UL2(i)=((1-rN(i)*the)*UL2(i)+ A(i))/(1+rN(i)*(1-the)); % Theta 1/2
     A=zeros(N,1); the=1/2;
      if(i>1) A(i)=A(i)+Neb1(i)*UL3(i-1); end
      if(i<N) A(i)=A(i)+Neb2(i)*UL3(i+1); end
     UL3(i)=((1-rN(i)*the)*UL3(i)+ A(i))/(1+rN(i)*(1-the)); % Theta
     A=zeros(N,1); the=1/2;
      if(i>1) A(i)=A(i)+Neb1(i)*UL4(i-1); end
      if(i<N) A(i)=A(i)+Neb2(i)*UL4(i+1); end
     UL4(i)=((1-rN(i)*the)*UL4(i)+ A(i))/(1+rN(i)*(1-the)); % Theta
      A=zeros(N,1); the=1/2;
      if(i>1) A(i)=A(i)+Neb1(i)*UL5(i-1); end
      if(i<N) A(i)=A(i)+Neb2(i)*UL5(i+1); end
     UL5(i)=((1-rN(i)*the)*UL5(i)+ A(i))/(1+rN(i)*(1-the)); % Theta
  end
    for i=1:2:N  %%% odd
       A=zeros(N,1);  
     if(i>1) A(i)=A(i)+Neb1(i)*UL1(i-1); end
     if(i<N) A(i)=A(i)+Neb2(i)*UL1(i+1); end
             if (tax(t+1)<tf-h)
    UL1(i)=UL1(i)*ee(i)+ A(i)/rN(i)*(1-ee(i));  % full
            else
    UL1(i)=UL1(i)*ee2(i)+ A(i)/rN(i)*(1-ee2(i)); %half 
             end
  A=zeros(N,1); the=1/2;
     if(i>1) A(i)=A(i)+Neb1(i)*UL2(i-1); end
     if(i<N) A(i)=A(i)+Neb2(i)*UL2(i+1); end
             if (tax(t+1)<tf-h)
     UL2(i)=((1-rN(i)*the)*UL2(i)+ A(i))/(1+rN(i)*(1-the)); % Theta 1/2 full
            else
     UL2(i)=((1-rN(i)/2*the)*UL2(i)+ A(i)/2)/(1+rN(i)/2*(1-the)); %half 
             end
   A=zeros(N,1); the=1/2;
     if(i>1) A(i)=A(i)+Neb1(i)*UL3(i-1); end
     if(i<N) A(i)=A(i)+Neb2(i)*UL3(i+1); end
             if (tax(t+1)<tf-h)
     UL3(i)=((1-rN(i)*the)*UL3(i)+ A(i))/(1+rN(i)*(1-the)); % Theta 1/2 full
            else
     UL3(i)=((1-rN(i)/2*the)*UL3(i)+ A(i)/2)/(1+rN(i)/2*(1-the)); %half 
             end
    A=zeros(N,1); the=1/2;
     if(i>1) A(i)=A(i)+Neb1(i)*UL4(i-1); end
     if(i<N) A(i)=A(i)+Neb2(i)*UL4(i+1); end
             if (tax(t+1)<tf-h)
     UL4(i)=((1-rN(i)*the)*UL4(i)+ A(i))/(1+rN(i)*(1-the)); % Theta 1/2 full
            else
     UL4(i)=((1-rN(i)/2*the)*UL4(i)+ A(i)/2)/(1+rN(i)/2*(1-the)); %half 
             end
    A=zeros(N,1); the=1/2;
     if(i>1) A(i)=A(i)+Neb1(i)*UL5(i-1); end
     if(i<N) A(i)=A(i)+Neb2(i)*UL5(i+1); end
             if (tax(t+1)<tf-h)
     UL5(i)=((1-rN(i)*the)*UL5(i)+ A(i))/(1+rN(i)*(1-the)); % Theta 1/2 full
            else
     UL5(i)=((1-rN(i)/2*the)*UL5(i)+ A(i)/2)/(1+rN(i)/2*(1-the)); %half 
            end
    end
 end
toc
%% L1 errors
Kul=zeros(N,1); MaxKul=0; 
for i=1:N
     Kul(i)=UL1(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,2)=MaxKul;
%% L2 errors
Kul=zeros(N,1); MaxKul=0; 
for i=1:N
     Kul(i)=UL2(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,3)=MaxKul;
%% L3 errors
Kul=zeros(N,1); MaxKul=0; 
for i=1:N
     Kul(i)=UL3(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,4)=MaxKul;%% L4 errors
Kul=zeros(N,1); MaxKul=0; 
for i=1:N
     Kul(i)=UL4(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,5)=MaxKul;
%% L5 errors
Kul=zeros(N,1); MaxKul=0; 
for i=1:N
     Kul(i)=UL5(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,6)=MaxKul;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=2*(T-1)+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF SIMULATION

Red	= [1, 0, 0];
L1c	 =    [123/255,	204/255,	181/255];
L2c	= [250/255,	50/255,		160/255]	;
L3c		= [0.5,	200/255,	255/255];

L4c =	[21/255,	27/255,		141/255];
L5c	=  [0.990, 0.8940, 0.250]; 
Heun = [200/255, 20/255, 10/255]; 

	%% figure 1:
  figure('Name', 'Errors as a fuction of time step h');   

	plot(Axhstep(:),MaxD(:,1), "--o", 'Color',Red, 'LineWidth', 2, 'MarkerSize',6); %% cn 
	hold on;
    plot(Axhstep(:),MaxD(:,2), '->',  'Color', L1c, 'LineWidth', 1.2, 'MarkerSize',5); %% cncn p=1/3
	hold on;
    plot(Axhstep(:),MaxD(:,3), '-p',  'Color', L2c, 'LineWidth', 1.2, 'MarkerSize',14); %% cncn p=2/3
	hold on;
    plot(Axhstep(:),MaxD(:,4), '-.s',  'Color', L3c, 'LineWidth', 2.3, 'MarkerSize',6); %% cncn p=1/2
	hold on;
	plot(Axhstep(:),MaxD(:,5), '-^',  'Color', L4c, 'LineWidth', 1.0,'MarkerSize',10); %% cncn p=3/4
	hold on;
    plot(Axhstep(:),MaxD(:,6), ':',  'Color', L5c, 'LineWidth', 4.6, 'MarkerSize',8); %% cncn p=3/4
	hold on;
      plot(Axhstep(:),MaxD(:,7), ':x',  'Color', Heun, 'LineWidth', 1.1, 'MarkerSize',22); %% cncn p=1/3
  	hold on;
	hold off;

XLabel = xlabel('time step size h');
set(XLabel, 'FontSize', 16);
set(XLabel, 'FontWeight', 'bold');

YLabel = ylabel('Errors');
set(YLabel, 'FontSize', 16); 
set(YLabel, 'FontWeight', 'bold');

Ax = gca;
Ax.XAxis.FontSize = 18;
Ax.YAxis.FontSize = 18;
Ax.FontWeight = 'bold';
Ax.TickLength = [0.018 0.035];
Ax.LineWidth = 1.2;

Legend=legend({'CNe','L1', 'L2','L3','L4','L5','Heun'},'Location','southeast', 'FontSize', 14);
set(Legend, 'FontSize', 14); 
set(Legend, 'FontWeight', 'bold');
set(gca,'xScale','log');
set(gca,'yScale','log');
grid on;



