%% Function File: Algorithm_1_v2_FullTheta.m 


	function Algorithm_1_v2_FullTheta( Label, Path_To_Result_Folder, Theta_0, Theta_1, Theta_2, Theta_3, Theta_4 )


		%%-------------------------------------------------------------------------------------------------------
		% STEP-0. Off a waning message.
			warning( 'off', 'MATLAB:xlswrite:AddSheet' ) ;

		%%-------------------------------------------------------------------------------------------------------
		% STEP-1. Display some information. 
			sRun = getRUN;
			NRun = int32( sRun.NRun );
			hThreshold = double( sRun.hThreshold  );
		
		
			Date_And_Time = datestr(now,'yyyy/mmmm/dd - HH:MM:SS.FFF ');
			message = sprintf("==> STRUC_1-A%i | Th0:%6.4f | Th1:%6.4f | Th2:%6.4f | Th3:%6.4f | Th4:%6.4f | START: %s \n", Label, Theta_0, Theta_1, Theta_2, Theta_3, Theta_4, Date_And_Time ); 
			fprintf( message );

		%%-------------------------------------------------------------------------------------------------------
		% STEP-2. Get data from the DTOA and ALGREF global variables. 
			sDtoA = getDTOA;
			AlgRef = getALGREF;

		%%-------------------------------------------------------------------------------------------------------
		% STEP-3. Set running parameters and algorithm parameters. 
		
			

		%%-------------------------------------------------------------------------------------------------------
		% STEP-4. Load data into local variables. 	
			
			Exist_Neighbour_Left   =  logical( sDtoA.Exist_Neighbour_Left );   
			Exist_Neighbour_Right  =  logical( sDtoA.Exist_Neighbour_Right );   
			Exist_Neighbour_Upper  =  logical( sDtoA.Exist_Neighbour_Upper );    
			Exist_Neighbour_Lower  =  logical( sDtoA.Exist_Neighbour_Lower );  
	
			Neighbour_Left  = double( sDtoA.Neighbour_Left ); 	
			Neighbour_Right = double( sDtoA.Neighbour_Right );	
			Neighbour_Upper = double( sDtoA.Neighbour_Upper );	
			Neighbour_Lower = double( sDtoA.Neighbour_Lower );	

			Nx    		= int32( sDtoA.Nx ); 
			N     		= int32( sDtoA.N ); 
			Mh    		= int32( sDtoA.Mh );
			
			timeInitial 	= double( sDtoA.timeInitial );
			timeFinal 		= double( sDtoA.timeFinal );
			timeInterval    = double( sDtoA.timeInterval ); 
			
			
			
			k     		= double( sDtoA.k );
			u0     		= double( sDtoA.u0 );	
			b          	= double( sDtoA.b );
			Axis_hstep 	= double( sDtoA.Axis_hstep );
		
			indexRow = int32( sDtoA.indexRow );
			indexColumn = int32( sDtoA.indexColumn );
			indexSum = int32( sDtoA.indexSum );
	
		%%-------------------------------------------------------------------------------------------------------
		% STEP-5. Definition of required local variables. 
		
			% --> Auxiliary variables.  
			a	     = double( 0 ); 
			uTemp 	 = zeros( N,   1,  'double' );  
			AbsDU    = zeros( N,   1,  'double' );  % AbsDU_i=|uREf_i-uNum_i| 
			% --> Storing results.
			uNum     = zeros( N,   1,  'double' );  % Temperatures obtained using the numerical algorithm at PhysT time.
			MaxD     = zeros( Mh,  1,  'double' );	% Maximum of the absolute values of the differences between the numerical and reference temperatures.	
			RunT 	 = zeros( Mh,  1,  'double' );  % The runtime of the algorithm in seconds at a given Time Step Size h.
			AveD     = zeros( Mh,  1,  'double' );
			DifQ     = zeros( Mh,  1,  'double' );  %	DifQ = KulQ   !!


		%%-------------------------------------------------------------------------------------------------------
		% STEP-6. Copmute the cell parities. We use the Hopscotch cell parity defition.

			cellParity =  int16( mod( indexSum, 2 ) ) ;
			for i=1:1:N
				if ( cellParity(i)==0 )
					cellParity(i) = 1; %%  EVEN cells. 
				else
					cellParity(i) = -1;%%  ODD cells. 
				end % if ( cellParity(i)==0 )
			end % for i=1:1:N

		
		%%-------------------------------------------------------------------------------------------------------
		% STEP-7. Implementation of the algorithm. 
				
			%% for loop over the Time Step Size variable h.
			for ih = 1:1:Mh
			
				h  = double( Axis_hstep(ih) );	
				e1 = double( exp(-b*h) ); 	


				hhalf = double( h/2 ); %Half time step size
				e1half = double( exp(-b*hhalf) );	
				
				if( h<hThreshold ) NRun = 1; end % if( h<hThreshold )  
				
					% Loop over the runs. Calculate the average of the RunT variable.
					for iRun=1:1:NRun 
						
						% Initialization of initial temperatures at time zero.
						uNum(:) = u0(:); 
						T = int32( ( timeInterval /h)+1 );	
						
						celltime = zeros( N, 1, 'double' );
						celltime(:) = timeInitial;
						
						%  MM = reshape(u0, Nx, Nz);
						%  disp( 'stop here' );
					
						% --> Start of Time Measurement. 
						tStart = tic;
						
						%% First h/2 time step for ODD cells.
						for i=1:1:N 
									
							if( cellParity(i)<0 ) % ODD cells			
									
								a = double( k(i) );
								celltime(i)=celltime(i)+hhalf;
						
								% Left Neighbour.
								if ( Exist_Neighbour_Left(i) )
									a = a + Neighbour_Left(i)*uNum(i-1); 
								end % if ( Exist_Neighbour_Left(i) )
									
								% Right Neighbour.
								if ( Exist_Neighbour_Right(i) )
									a = a + Neighbour_Right(i)*uNum(i+1);  
								end % if ( Exist_Neighbour_Right(i) )
							
								% Upper Neighbour.
								if ( Exist_Neighbour_Upper(i) )
									a = a + Neighbour_Upper(i)*uNum(i-Nx);		
								end % if ( Exist_Neighbour_Upper(i) )
							
								% Lower Neighbour.
								if ( Exist_Neighbour_Lower(i) )
									a = a + Neighbour_Lower(i)*uNum(i+Nx);
								end % if ( Exist_Neighbour_Lower(i) )
							
								% 0 - MODIFY THIS LINE! 	
								uNum(i) = ((1-b(i)*hhalf*Theta_0)*uNum(i) + hhalf*a)/(1+(b(i)*hhalf*(1-Theta_0)));
	
							end % if( cellParity(i)>0 ) ODD cells	
								
						end % for i=1:1:N
	
                        
						% for loop over the Physical Time.
						for t = 1:2:(T-2)		
                            
							% Loop over the cells.
							for i = 1:1:N
								
								if( cellParity(i)>0 ) % EVEN cells
								
									a = double( k(i) ); 
									celltime(i)=celltime(i)+h;
							
									% Left Neighbour.
									if ( Exist_Neighbour_Left(i) )
										a = a + Neighbour_Left(i)*uNum(i-1); 
									end % if ( Exist_Neighbour_Left(i) )
										
									% Right Neighbour.
									if ( Exist_Neighbour_Right(i) )
										a = a + Neighbour_Right(i)*uNum(i+1);  
									end % if ( Exist_Neighbour_Right(i) )
								
									% Upper Neighbour.
									if ( Exist_Neighbour_Upper(i) )
										a = a + Neighbour_Upper(i)*uNum(i-Nx);		
									end % if ( Exist_Neighbour_Upper(i) )
								
									% Lower Neighbour.
									if ( Exist_Neighbour_Lower(i) )
										a = a + Neighbour_Lower(i)*uNum(i+Nx);
									end % if ( Exist_Neighbour_Lower(i) )
										
									uNum(i) = ((1-b(i)*h*Theta_1)*uNum(i) + h*a)/(1+(b(i)*h*(1-Theta_1)));								
									
								end % if( cellParity(i)<0 ) % EVEN cells
									
							end % for i = 1:1:N

							% Loop over the cells.
							for i = 1:1:N
								
								if(cellParity(i)<0) % ODD cells 
								
									a = double( k(i) );
									celltime(i)=celltime(i)+h;
							
									% Left Neighbour.
									if ( Exist_Neighbour_Left(i) )
										a = a + Neighbour_Left(i)*uNum(i-1); 
									end % if ( Exist_Neighbour_Left(i) )
										
									% Right Neighbour.
									if ( Exist_Neighbour_Right(i) )
										a = a + Neighbour_Right(i)*uNum(i+1);  
									end % if ( Exist_Neighbour_Right(i) )
								
									% Upper Neighbour.
									if ( Exist_Neighbour_Upper(i) )
										a = a + Neighbour_Upper(i)*uNum(i-Nx);		
									end % if ( Exist_Neighbour_Upper(i) )
								
									% Lower Neighbour.
									if ( Exist_Neighbour_Lower(i) )
										a = a + Neighbour_Lower(i)*uNum(i+Nx);
									end % if ( Exist_Neighbour_Lower(i) )
										
									uNum(i) = ((1-b(i)*h*Theta_2)*uNum(i) + h*a)/(1+(b(i)*h*(1-Theta_2)));
											
								end % if( cellParity(i)>0 ) ODD cells
									
							end % for i = 1:1:N


							% Loop over the cells.
							for i = 1:1:N
								
								if(cellParity(i)>0) % EVEN cells
								
									a = double( k(i) );
									celltime(i)=celltime(i)+h;
							
									% Left Neighbour.
									if ( Exist_Neighbour_Left(i) )
										a = a + Neighbour_Left(i)*uNum(i-1); 
									end % if ( Exist_Neighbour_Left(i) )
										
									% Right Neighbour.
									if ( Exist_Neighbour_Right(i) )
										a = a + Neighbour_Right(i)*uNum(i+1);  
									end % if ( Exist_Neighbour_Right(i) )
								
									% Upper Neighbour.
									if ( Exist_Neighbour_Upper(i) )
										a = a + Neighbour_Upper(i)*uNum(i-Nx);		
									end % if ( Exist_Neighbour_Upper(i) )
								
									% Lower Neighbour.
										if ( Exist_Neighbour_Lower(i) )
										a = a + Neighbour_Lower(i)*uNum(i+Nx);
									end % if ( Exist_Neighbour_Lower(i) )
	
									uNum(i) = ((1-b(i)*h*Theta_3)*uNum(i) + h*a)/(1+(b(i)*h*(1-Theta_3)));
											
								end % if( cellParity(i)>0 ) EVEN cells
									
							end % for i = 1:1:N	


							% Loop over the cells.
							for i = 1:1:N
								
								if( cellParity(i)<0 ) % ODD cells
								
									a = double( k(i) );
									celltime(i)=celltime(i)+h;
									
									% Left Neighbour.
									if ( Exist_Neighbour_Left(i) )
										a = a + Neighbour_Left(i)*uNum(i-1); 
									end % if ( Exist_Neighbour_Left(i) )
										
									% Right Neighbour.
									if ( Exist_Neighbour_Right(i) )
										a = a + Neighbour_Right(i)*uNum(i+1);  
									end % if ( Exist_Neighbour_Right(i) )
								
									% Upper Neighbour.
									if ( Exist_Neighbour_Upper(i) )
										a = a + Neighbour_Upper(i)*uNum(i-Nx);		
									end % if ( Exist_Neighbour_Upper(i) )
								
									% Lower Neighbour.
									if ( Exist_Neighbour_Lower(i) )
										a = a + Neighbour_Lower(i)*uNum(i+Nx);
									end % if ( Exist_Neighbour_Lower(i) )
									
									% 4 - MODIFY THIS LINE!
									if ( celltime(i)<timeInterval ) % Make an h step.
											uNum(i) = ((1-b(i)*h*Theta_4)*uNum(i) + h*a)/(1+(b(i)*h*(1-Theta_4)));
									else % Make an hhalf step.  
										uNum(i) = ((1-b(i)*hhalf*Theta_4)*uNum(i) + hhalf*a)/(1+(b(i)*hhalf*(1-Theta_4)));
										
										celltime(i)=celltime(i)-hhalf;
									end % if ( celltime(i)<PhysT )

								end % if( cellParity(i)<0 ) % ODD cells
								
							end % for i = 1:1:N
							
							
						end % for t = 1:2:(T-2)

						% --> Stop of Time Measurement. 
						tStop = toc( tStart );
						% Evaluate Running Time.
						RunT(ih) = RunT(ih)+tStop/double( NRun ); 
						
							
					end % for iRun=1:1:NRun 


				%------------------------------------------------------------------------------------------------------
				% Compute errors for each Time Step Size.
				% Compute AbsDU.
				
				% Reference: uAN !! 
				AbsDU = abs(sDtoA.uRef - uNum );
				
				% Compute MaxD. 
				MaxD(ih) = max ( AbsDU );
				% Compute AveD.
				AveD(ih) = mean( AbsDU );
				% Compute DifQ.
				DifQ(ih) = sum( sDtoA.C .* AbsDU );
			% end % if( hThreshold )
			
		end % for ih=1:1:Mh - Loop over the time step size h. 

		%%-------------------------------------------------------------------------------------------------------
		% STEP-8. Calculate some aggregated error.
		MaxD_REF = double( AlgRef.MaxD_REF ); 
		AveD_REF = double( AlgRef.AveD_REF ); 
		DifQ_REF = double( AlgRef.DifQ_REF );
		
		Sum_Log_MaxD = double( sum(log10( MaxD ))); Sum_Log_MaxD_REF = double( sum(log10( MaxD_REF )));		
		Sum_Log_AveD = double( sum(log10( AveD ))); Sum_Log_AveD_REF = double( sum(log10( AveD_REF )));
		Sum_Log_DifQ = double( sum(log10( DifQ ))); Sum_Log_DifQ_REF = double( sum(log10( DifQ_REF )));

		Rel_Err_MaxD = double( Sum_Log_MaxD / Sum_Log_MaxD_REF ); 
		Rel_Err_AveD = double( Sum_Log_AveD / Sum_Log_AveD_REF ); 
		Rel_Err_DifQ = double( Sum_Log_DifQ / Sum_Log_DifQ_REF ); 

		Delta_Err_MaxD = double( Sum_Log_MaxD_REF - Sum_Log_MaxD )/ double( Mh ); 
		Delta_Err_AveD = double( Sum_Log_AveD_REF - Sum_Log_AveD )/ double( Mh );
		Delta_Err_DifQ = double( Sum_Log_DifQ_REF - Sum_Log_DifQ )/ double( Mh );
		
		Av_Error = double( ( Delta_Err_MaxD + Delta_Err_AveD + Delta_Err_DifQ ) / ( 3.0 ) );  


		%%-------------------------------------------------------------------------------------------------------
		% STEP-9. Write data to output file.  
		
			% Part-1.
				Current_File = fullfile(Path_To_Result_Folder, "Algorithm_Table_All.xls");	
				Wrote_Sheet = sprintf('A%i', Label);

				Header = ["h" "MaxD" "AveD" "DifQ" "RunT" ];	
				writematrix(Header , Current_File, 'Sheet', Wrote_Sheet, 'Range', 'A1:E1' );	
				
				Wrote_Range = sprintf('A2:A%i',Mh+1);
				writematrix(Axis_hstep, Current_File, 'Sheet', Wrote_Sheet, 'Range', Wrote_Range ); 
			
				Wrote_Range = sprintf('B2:B%i',Mh+1);
				writematrix(MaxD, Current_File, 'Sheet', Wrote_Sheet, 'Range', Wrote_Range ); 

				Wrote_Range = sprintf('C2:C%i',Mh+1);
				writematrix(AveD, Current_File, 'Sheet', Wrote_Sheet, 'Range', Wrote_Range ); 

				Wrote_Range = sprintf('D2:D%i',Mh+1);
				writematrix(DifQ, Current_File, 'Sheet', Wrote_Sheet, 'Range', Wrote_Range ); 

				Wrote_Range = sprintf('E2:E%i',Mh+1);
				writematrix(RunT, Current_File, 'Sheet', Wrote_Sheet, 'Range', Wrote_Range ); 

			% Part-2. 
				if ( Label == 1 )
				
					% Table_MaxD
					Current_File = fullfile(Path_To_Result_Folder, "Algorithm_Table_MaxD.xls");					
					Wrote_Sheet = "MaxD";
					writematrix("h", Current_File, 'Sheet', Wrote_Sheet, 'Range', "A1:A1" );
					Wrote_Range = sprintf('A2:A%i',Mh+1);
					writematrix(Axis_hstep, Current_File, 'Sheet', Wrote_Sheet, 'Range', Wrote_Range );
				
					% Table_AveD
					Current_File = fullfile(Path_To_Result_Folder, "Algorithm_Table_AveD.xls");					
					Wrote_Sheet = "AveD";
					writematrix("h", Current_File, 'Sheet', Wrote_Sheet, 'Range', "A1:A1" );
					Wrote_Range = sprintf('A2:A%i',Mh+1);
					writematrix(Axis_hstep, Current_File, 'Sheet', Wrote_Sheet, 'Range', Wrote_Range );
				
					% Table_DifQ
					Current_File = fullfile(Path_To_Result_Folder, "Algorithm_Table_DifQ.xls");					
					Wrote_Sheet = "DifQ";
					writematrix("h", Current_File, 'Sheet', Wrote_Sheet, 'Range', "A1:A1" );
					Wrote_Range = sprintf('A2:A%i',Mh+1);
					writematrix(Axis_hstep, Current_File, 'Sheet', Wrote_Sheet, 'Range', Wrote_Range );
				
					% Table_DifQ
					Current_File = fullfile(Path_To_Result_Folder, "Algorithm_Table_RunT.xls");					
					Wrote_Sheet = "RunT";
					writematrix("h", Current_File, 'Sheet', Wrote_Sheet, 'Range', "A1:A1" );
					Wrote_Range = sprintf('A2:A%i',Mh+1);
					writematrix(Axis_hstep, Current_File, 'Sheet', Wrote_Sheet, 'Range', Wrote_Range );						
				
				end % if ( Label == 1 ) 
				
				% Get the right letter. 
				Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
				Letter = Alphabet( Label+1 );

				% Table_MaxD.
				Current_File = fullfile(Path_To_Result_Folder, "Algorithm_Table_MaxD.xls");					
				Wrote_Sheet = "MaxD";
				Wrote_Range = sprintf('%s1:%s1',Letter, Letter);
				Header = sprintf("A%i", Label);
				writematrix(Header, Current_File, 'Sheet', Wrote_Sheet, 'Range', Wrote_Range );
				Wrote_Range = sprintf('%s2:%s%i',Letter, Letter, Mh+1);
				writematrix(MaxD, Current_File, 'Sheet', Wrote_Sheet, 'Range', Wrote_Range );

				% Table_AveD.
				Current_File = fullfile(Path_To_Result_Folder, "Algorithm_Table_AveD.xls");					
				Wrote_Sheet = "AveD";
				Wrote_Range = sprintf('%s1:%s1',Letter, Letter);
				Header = sprintf("A%i", Label);
				writematrix(Header, Current_File, 'Sheet', Wrote_Sheet, 'Range', Wrote_Range );
				Wrote_Range = sprintf('%s2:%s%i',Letter, Letter, Mh+1);
				writematrix(AveD, Current_File, 'Sheet', Wrote_Sheet, 'Range', Wrote_Range );

				% Table_DifQ.
				Current_File = fullfile(Path_To_Result_Folder, "Algorithm_Table_DifQ.xls");					
				Wrote_Sheet = "DifQ";
				Wrote_Range = sprintf('%s1:%s1',Letter, Letter);
				Header = sprintf("A%i", Label);
				writematrix(Header, Current_File, 'Sheet', Wrote_Sheet, 'Range', Wrote_Range );
				Wrote_Range = sprintf('%s2:%s%i',Letter, Letter, Mh+1);
				writematrix(DifQ, Current_File, 'Sheet', Wrote_Sheet, 'Range', Wrote_Range );
				
				% Table_RunT.
				Current_File = fullfile(Path_To_Result_Folder, "Algorithm_Table_RunT.xls");					
				Wrote_Sheet = "RunT";
				Wrote_Range = sprintf('%s1:%s1',Letter, Letter);
				Header = sprintf("A%i", Label);
				writematrix(Header, Current_File, 'Sheet', Wrote_Sheet, 'Range', Wrote_Range );
				Wrote_Range = sprintf('%s2:%s%i',Letter, Letter, Mh+1);
				writematrix(RunT, Current_File, 'Sheet', Wrote_Sheet, 'Range', Wrote_Range );

			% Part-3.
				if ( Label == 1 )
					Current_File = fullfile(Path_To_Result_Folder, "Algorithm_Table_Temp_hMin.xls");	
					Wrote_Sheet = "Temperatures";
					
					writematrix("Cell index", Current_File, 'Sheet', Wrote_Sheet, 'Range', "A1:A1" );
					Wrote_Range = sprintf('A2:A%i',N+1);
					i = transpose( int32( 1:1:N ) );
					writematrix(i, Current_File, 'Sheet', Wrote_Sheet, 'Range', Wrote_Range );
					
				end % if ( Label == 1 )

				% Table_uNum.
				Current_File = fullfile(Path_To_Result_Folder, "Algorithm_Table_Temp_hMin.xls");					
				Wrote_Sheet = "Temperatures";
				Wrote_Range = sprintf('%s1:%s1',Letter, Letter);
				Header = sprintf("A%i", Label);
				writematrix(Header, Current_File, 'Sheet', Wrote_Sheet, 'Range', Wrote_Range );
				Wrote_Range = sprintf('%s2:%s%i',Letter, Letter, N+1);
				writematrix(uNum, Current_File, 'Sheet', Wrote_Sheet, 'Range', Wrote_Range );

			% Part-4.
				Current_File = fullfile(Path_To_Result_Folder, "Algorithm_List.xls");	
				Wrote_Range = sprintf('A%i:A%i', Label, Label );
				writematrix(message, Current_File, 'Sheet', 'List', 'Range', Wrote_Range );


			% Part-5.  
				Current_File = fullfile(Path_To_Result_Folder, "Algorithm_Table_Aggregated.xls");
				Wrote_Sheet = "AggrErrors";
				if( Label == 1 )
					Header = ["Algorithm"  "Av_Error" "Delta_Err_MaxD" "Delta_Err_AveD" "Delta_Err_DifQ" "Rel_Err_MaxD" "Rel_Err_AveD" "Rel_Err_DifQ"  "Sum_Log_MaxD" "Sum_Log_AveD" "Sum_Log_DifQ" ];
					Wrote_Range = "A1:J1";
					writematrix(Header , Current_File, 'Sheet', Wrote_Sheet, 'Range', Wrote_Range );	
				end % if( Label==1 )
				
				message = sprintf(" STRUC1-A%i | Th0:%6.4f | Th1:%6.4f | Th2:%6.4f | Th3:%6.4f | Th4:%6.4f | \n", Label, Theta_0, Theta_1, Theta_2, Theta_3, Theta_4); 
				Saved_Data = [message Av_Error Delta_Err_MaxD Delta_Err_AveD Delta_Err_DifQ Rel_Err_MaxD Rel_Err_AveD Rel_Err_DifQ Sum_Log_MaxD Sum_Log_AveD Sum_Log_DifQ ];
				Wrote_Range = sprintf("A%i:J%i", Label+1, Label+1);
				writematrix(Saved_Data , Current_File, 'Sheet', Wrote_Sheet, 'Range', Wrote_Range );
				
                % disp('itt');
			
	end % Algorithm_1_v2_FullTheta


