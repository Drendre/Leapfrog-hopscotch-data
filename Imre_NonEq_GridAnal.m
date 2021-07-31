%% Test against Imre's exact solution in 1D, 
%% Generating non-equidistant grid   
clear;close all;fprintf('\n')
D=1; %How many different timestep we want 
%%%%%%%  Generating Grid
dxZero = double( 0.01 ); % 0.015 
gamma = double( 1e-11 ); % 3e-6
Nx_Coord_Pos = int32( 1000 ); % 500
x_Coord_Pos = zeros( Nx_Coord_Pos, 1, 'double' );
 x_Coord_Pos(1) = dxZero; 
	for i=2:1:Nx_Coord_Pos
	x_Coord_Pos(i) = x_Coord_Pos(i-1) + exp( gamma*( double(i)-1 )^4 )*dxZero; %^4
	end % for i=1:1:Nx	
 	xCoord = vertcat([-flip(x_Coord_Pos); 0], x_Coord_Pos); 
 	Nx_Coord = numel( xCoord ); 
 	N = int32( Nx_Coord-1 );
 	dx = zeros( N, 1, 'double' ); 
Cent = zeros( N, 1, 'double' ); %cell centers
R = zeros( N-1, 1, 'double' );
C = zeros( N, 1, 'double' );
		for i=1:1:N
			dx(i) = abs( xCoord(i+1)-xCoord(i) ); %lenght of cells
            C(i)=dx(i);
            Cent(i) = xCoord(i)+ dx(i)/2;
       	end % for i=1:1:Nx
        for i=1:1:N-1
	        R(i) = abs( Cent(i+1)-Cent(i) ); %distance of cell-centers
 		end % for i=1:1:Nx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ti=0.5; tf=1; TIME=tf-ti; %Diffusion coefficient, lenght of the space and time interval, x parameter
M1=zeros(N,N); M1(1,1)=-1/R(1)/C(1);M1(1,2)=1/R(1)/C(1);M1(N,N)=-1/R(N-1)/C(N);M1(N,N-1)=1/R(N-1)/C(N);
for i=2:N-1 
    M1(i,i-1)=1/R(i-1)/C(i); M1(i,i)=-1*(1/R(i)+1/R(i-1))/C(i); M1(i,i+1)=1/R(i)/C(i);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u0=zeros(N,1); uex=zeros(N,1); %%initial values and exact solution
for i=1:N 
 u0(i)=Cent(i)/ti^(5/2)*exp(-Cent(i)^2/4/D/ti)*(1-Cent(i)^2/6/D/ti); 
 uex(i)=Cent(i)/tf^(5/2)*exp(-Cent(i)^2/4/D/tf)*(1-Cent(i)^2/6/D/tf);  %%Initial and exact final results
end
plot(Cent,u0,Cent,uex,':')