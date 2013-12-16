//--------------------------------------------------------------------------------



// Disabled/null histogram_option for prescan
if(histogram_option==0){ 
}



//--------------------------------------------------------------------------------



// sort 1d spectra of parameters relevant to Gammasphere and the data i.e. not Microball
if(histogram_option==1){ 

   for(i=1; i<=clean_ge; i++){
      ge_id[i]=fmod(ge_id[i],1000);
      if(ge_energy[i]/gscomp[0] >= 0 && ge_energy[i]/gscomp[0] < SPEC_LENGTH){
         spec[ge_id[i]][(short)(ge_energy[i]/gscomp[0])]++;
      }
      if(ge_time[i]/gscomp[2] >= 0 && ge_time[i]/gscomp[2] < SPEC_LENGTH){
         spec[ge_id[i]+120][(short)(ge_time[i]/gscomp[2])]++;
      }
      spec[240][ge_id[i]]++;
      spec[250][(short)(ge_energy[i]/gscomp[0])]++;
   }
				
   if(clean_ge && clean_ge<SPEC_LENGTH) 			spec[241][clean_ge]++;
   if(dirty_ge && dirty_ge<SPEC_LENGTH) 			spec[242][dirty_ge]++;
   if(bgo_only && bgo_only<SPEC_LENGTH) 			spec[243][bgo_only]++;
   if(ge_low_res_sum >= 0 && ge_low_res_sum < SPEC_LENGTH) 	spec[244][ge_low_res_sum]++;
   if(bgo_low_res_sum >= 0 && bgo_low_res_sum < SPEC_LENGTH) 	spec[245][bgo_low_res_sum]++;
   if(rf >= 0 && rf < SPEC_LENGTH) 				spec[246][rf]++;
   if(total_fold && total_fold < SPEC_LENGTH) 			spec[247][total_fold]++;
   if(tac1>=0 && tac1 < SPEC_LENGTH) 				spec[248][tac1]++;
   if(event_length && event_length<SPEC_LENGTH) 		spec[249][event_length]++;

}



//--------------------------------------------------------------------------------



// Stripped down version of histogram option 1 with no time spectra
// Set up to sort source spectra at 1/3 keV per channel and accommodate the fact that 
// gcc apparently can't cope with 300 8192 channel spectra
if(histogram_option==101){	

   for(i=1; i<=clean_ge; i++){
      ge_id[i]=fmod(ge_id[i],1000);
		
      if(ge_energy[i]/gscomp[0] >= 0 && ge_energy[i]/gscomp[0] < SPEC_LENGTH){
         spec[ge_id[i]][(short)(ge_energy[i]/gscomp[0])]++;
      }
      spec[140][ge_id[i]]++;
      spec[160][(short)(ge_energy[i]/gscomp[0])]++;
   }
				
   if(clean_ge && clean_ge<SPEC_LENGTH)
      spec[141][clean_ge]++;
   if(dirty_ge && dirty_ge<SPEC_LENGTH)
      spec[142][dirty_ge]++;
   if(bgo_only && bgo_only<SPEC_LENGTH)
      spec[143][bgo_only]++;
   if(ge_low_res_sum >= 0 && ge_low_res_sum < SPEC_LENGTH)
      spec[144][ge_low_res_sum]++;
   if(bgo_low_res_sum >= 0 && bgo_low_res_sum < SPEC_LENGTH)
      spec[145][bgo_low_res_sum]++;
   if(rf >= 0 && rf < SPEC_LENGTH)
      spec[146][rf]++;
   if(tac1>=0 && tac1 < SPEC_LENGTH)
      spec[148][tac1]++;
   if(total_fold && total_fold < SPEC_LENGTH)
      spec[147][total_fold]++;
   if(ppac_x > 0 && ppac_x < SPEC_LENGTH)
      spec[150][ppac_x]++;
   if(clean_ge==0 && ppac_x > 0 && ppac_x < SPEC_LENGTH)
      spec[151][ppac_x]++;
   if(clean_ge>0 && ppac_x > 0 && ppac_x < SPEC_LENGTH)
      spec[152][ppac_x]++;
   if(event_length && event_length<SPEC_LENGTH)
      spec[153][event_length]++;

}



//--------------------------------------------------------------------------------


// Sort ring-by-ring 1D spectra

if(histogram_option==102){ 

   for(i=1; i<=clean_ge; i++){
      ge_id[i]=fmod(ge_id[i],1000);
      if(ge_energy[i]/gscomp[0] >= 0 && ge_energy[i]/gscomp[0] < SPEC_LENGTH){
         spec[gs_ring[ge_id[i]]][(short)(ge_energy[i]/gscomp[0])]++;
      }
      if(ge_time[i]/gscomp[2] >= 0 && ge_time[i]/gscomp[2] < SPEC_LENGTH){
         spec[gs_ring[ge_id[i]+120]][(short)(ge_time[i]/gscomp[2])]++;
      }
      spec[240][ge_id[i]]++;
      spec[250][(short)(ge_energy[i]/gscomp[0])]++;
   }
				
   if(clean_ge && clean_ge<SPEC_LENGTH) 			spec[241][clean_ge]++;
   if(dirty_ge && dirty_ge<SPEC_LENGTH)  			spec[242][dirty_ge]++;
   if(bgo_only && bgo_only<SPEC_LENGTH) 			spec[243][bgo_only]++;
   if(ge_low_res_sum >= 0 && ge_low_res_sum < SPEC_LENGTH) 	spec[244][ge_low_res_sum]++;
   if(bgo_low_res_sum >= 0 && bgo_low_res_sum < SPEC_LENGTH) 	spec[245][bgo_low_res_sum]++;
   if(rf >= 0 && rf < SPEC_LENGTH) 				spec[246][rf]++;
   if(tac1>=0 && tac1 < SPEC_LENGTH) 				spec[247][tac1]++;
   if(total_fold && total_fold < SPEC_LENGTH) 			spec[248][total_fold]++;
   if(event_length && event_length<SPEC_LENGTH) 		spec[249][event_length]++;

}
      


//------------------------------------------------------------------------------



// clean Ge side channel information
if(histogram_option==103){
			
   for(i=1;i<=clean_ge;i++){
      ge_id[i]=fmod(ge_id[i],1000);
      if(ge_side[i] && ge_side[i] <SPEC_LENGTH) spec[ge_id[i]][ge_side[i]]++;
   }
}



//--------------------------------------------------------------------------------



// sort Gammasphere energy against usec clock maps
if(histogram_option==2){

   for(i=1;i<=clean_ge;i++){ 
      ge_id[i]=fmod(ge_id[i],1000);
      if(ge_energy[i]/gscomp[0] < MAT_LENGTH-1 && ge_energy[i]/gscomp[0] > 0){
         matrix[ge_id[i]][(int)(usec_clock/500e6+(float)rand()/RMAX-0.5)][(short)(ge_energy[i]/gscomp[0]+(float)rand()/RMAX-0.5)]++;
      }
   }
}



//------------------------------------------------------------------------------



// sort Gammasphere time against usec clock maps
if(histogram_option==201){

   for(i=1;i<=clean_ge;i++){ 
      ge_id[i]=fmod(ge_id[i],1000);
      if(ge_time[i]/gscomp[2] < MAT_LENGTH-1 && ge_time[i]/gscomp[2] > 0){
         matrix[ge_id[i]][(int)(usec_clock/500e6+(float)rand()/RMAX-0.5)][(short)(ge_time[i]/gscomp[2]+(float)rand()/RMAX-0.5)]++;
      }
   }
}



//------------------------------------------------------------------------------


// sort Microball energy against usec clock
if(histogram_option==202 && microball){

   for(i=1;i<=mball_rank_fold[2];i++){ 
      /* Remember: 	
         mb_pass[i][1]=energy
         mb_pass[i][2]=time
         mb_pass[i][3]=pid
         mb_pass[i][4]=id  
         mb_pass[i][5]=ratio */
      if(mb_pass[i][1]/mbcomp[0]<MAT_LENGTH && mb_pass[i][1]/mbcomp[0]>0){
         matrix[mb_pass[i][4]][(int)(usec_clock/500e6+(float)rand()/RMAX-0.5)][(short)(mb_pass[i][1]/mbcomp[0]+(float)rand()/RMAX-0.5)]++;
      }
   }
}



//------------------------------------------------------------------------------


// sort Microball time against usec clock
if(histogram_option==203 && microball){

   for(i=1;i<=mball_rank_fold[2];i++){ 
      /* Remember: 	
         mb_pass[i][1]=energy
         mb_pass[i][2]=time
         mb_pass[i][3]=pid
         mb_pass[i][4]=id  
         mb_pass[i][5]=ratio */

      if(mb_pass[i][2]/mbcomp[2]<MAT_LENGTH && mb_pass[i][2]/mbcomp[2]>0){
         matrix[mb_pass[i][4]][(int)(usec_clock/500e6+(float)rand()/RMAX-0.5)][(short)(mb_pass[i][2]/mbcomp[2]+(float)rand()/RMAX-0.5)]++;
      }
   }
}



//--------------------------------------------------------------------------------



// sort file # against usec clock maps
if(histogram_option==204){

   for(i=1;i<=clean_ge;i++){ 
      ge_id[i]=fmod(ge_id[i],1000);
      if(usec_clock/500e6<MAT_LENGTH && usec_clock>0){
         matrix[0][(int)(usec_clock/500e6+(float)rand()/RMAX-0.5)][counter[3]]++;
      }
   }
}



//------------------------------------------------------------------------------



// sort Microball energy, time, PID vs ID spectra
if(histogram_option==3 && microball){

   //printf("Event number: %d\n",counter[8]);
   for(i=1;i<=mball_rank_fold[2];i++){ 
      /* Remember mball_rank_fold[2] is the mball fold when there is a PID, 
         energy and time all present for all CsI detectors */

      /* Remember: 	
         mb_pass[i][1]=energy
         mb_pass[i][2]=time
         mb_pass[i][3]=pid
         mb_pass[i][4]=id  
         mb_pass[i][5]=ratio */

      spec[mb_pass[i][4]    ][(short)(mb_pass[i][1]/mbcomp[0]+(float)rand()/RMAX-0.5)]++;
      spec[mb_pass[i][4]+100][(short)(mb_pass[i][2]/mbcomp[2]+(float)rand()/RMAX-0.5)]++;
      spec[mb_pass[i][4]+200][(short)(mb_pass[i][3]/mbcomp[1]+(float)rand()/RMAX-0.5)]++;
      spec[300][mb_pass[i][4]]++;
   }
	
   spec[301][mball_rank_fold[0]]++;
   spec[302][mball_rank_fold[1]]++;
   spec[303][mball_rank_fold[2]]++;

}



//------------------------------------------------------------------------------




// sort run-by-run Microball spectra
if(histogram_option==301 && microball){

   for(i=1;i<=mball_rank_fold[2];i++){ 
      /* Remember: 	
         mb_pass[i][1]=energy
         mb_pass[i][2]=time
         mb_pass[i][3]=pid
         mb_pass[i][4]=id  
         mb_pass[i][5]=ratio */
      if(counter[3]<SPEC_LENGTH && mb_pass[i][2]/mbcomp[2]<SPEC_LENGTH && mb_pass[i][2]>0){
         spec[counter[3]][(short)(mb_pass[i][2]/mbcomp[2]+(float)rand()/RMAX-0.5)]++;
      }
   }
}



//------------------------------------------------------------------------------


// sort Microball energy vs PID maps
if(histogram_option==4 && microball){ 

   for(i=1;i<=mball_rank_fold[2];i++){
      if(mb_pass[i][1]/mbcomp[0]<MAT_LENGTH && mb_pass[i][3]/mbcomp[1]<MAT_LENGTH && mb_pass[i][1]/mbcomp[0]>0 && mb_pass[i][3]/mbcomp[1]>0){
         //printf("mbcomp[0] %f mbcomp[1] %f mbcomp[2] %f\n",mbcomp[0],mbcomp[1],mbcomp[2]);
         //printf("id %d x %d y %d\n",mb_pass[i][4],(short)(mb_pass[i][1]/mbcomp[0]),(short)(mb_pass[i][3]/mbcomp[1]));
         //spec[1][mb_pass[i][4]]++;
         matrix[mb_pass[i][4]][(short)(mb_pass[i][1]/mbcomp[0]+(float)rand()/RMAX-0.5)][(short)(mb_pass[i][3]/mbcomp[1]+(float)rand()/RMAX-0.5)]++;
      }
   }
}



//------------------------------------------------------------------------------



// sort Microball energy vs time maps
if(histogram_option==5){ 

   for(i=1;i<=mball_rank_fold[2];i++){
      if(mb_pass[i][2]/mbcomp[2] < MAT_LENGTH && mb_pass[i][1]/mbcomp[0] < MAT_LENGTH){
         matrix[mb_pass[i][4]][(short)(mb_pass[i][1]/mbcomp[0]+(float)rand()/RMAX-0.5)][(short)(mb_pass[i][2]/mbcomp[2]+(float)rand()/RMAX-0.5)]++;
      }
   }
}



//------------------------------------------------------------------------------




// sort Microball energy vs ratio maps
if(histogram_option==6){ 

   for(i=1;i<=mball_rank_fold[2];i++){
      if(mb_pass[i][1]/mbcomp[0] > 0 && mb_pass[i][1]/mbcomp[0] < MAT_LENGTH && mb_pass[i][5] > 0 && mb_pass[i][5] < MAT_LENGTH){
         matrix[mb_pass[i][4]][(short)(mb_pass[i][5]+(float)rand()/RMAX-0.5)][(short)(mb_pass[i][1]/mbcomp[0]+(float)rand()/RMAX-0.5)]++;
      }
   }
}



//------------------------------------------------------------------------------



// sort time vs ratio maps
if(histogram_option==7){ 

   for(i=1;i<=mball_rank_fold[2];i++){
      if(mb_pass[i][2]/mbcomp[2] > 0 && mb_pass[i][2]/mbcomp[2] < MAT_LENGTH && mb_pass[i][5] > 0 && mb_pass[i][5] < MAT_LENGTH){
         matrix[mb_pass[i][4]][(short)(mb_pass[i][5]+(float)rand()/RMAX-0.5)][(short)(mb_pass[i][2]/mbcomp[2]+(float)rand()/RMAX-0.5)]++;
      }
   }
}

//------------------------------------------------------------------------------



// sort a symmetrized ungated gamma-gamma matrix with RMAX
if(histogram_option==8){ 

   for(i=1; i<=clean_ge; i++){
      if(ge_energy[i]/gscomp[0]>10 && ge_energy[i]/gscomp[0] < MAT_LENGTH-1){
         for(j=1; j<=clean_ge; j++){
            if(i!=j && ge_energy[j]/gscomp[0]>10 && ge_energy[j]/gscomp[0] < MAT_LENGTH-1){
               matrix[0][(short)(ge_energy[i]/gscomp[0]+(float)rand()/RMAX-0.5)][(short)(ge_energy[j]/gscomp[0]+(float)rand()/RMAX-0.5)]++;
            }
         }
      }
   }
}



//------------------------------------------------------------------------------



// sort a symmetrized ungated gamma gamma matrix with RAND_MAX
if(histogram_option==801){ 

   for(i=1; i<=clean_ge; i++){
      if(ge_energy[i]/gscomp[0]>10 && ge_energy[i]/gscomp[0] < MAT_LENGTH-1){
         for(j=1; j<=clean_ge; j++){
            if(i!=j && ge_energy[j]/gscomp[0]>0 && ge_energy[j]/gscomp[0] < MAT_LENGTH-1){
               matrix[0][(short)(ge_energy[i]/gscomp[0]+(float)rand()/RAND_MAX-0.5)][(short)(ge_energy[j]/gscomp[0]+(float)rand()/RAND_MAX-0.5)]++;
            }
         }
      }
   }
}





//------------------------------------------------------------------------------



// sort a Gammasphere energy-time matrix
if(histogram_option==802){
  for(i=1; i<=clean_ge; i++){
    if(ge_energy[i]/gscomp[0] < MAT_LENGTH-1 && ge_time[i]/gscomp[2] < MAT_LENGTH-1){
      matrix[0][(short)(ge_energy[i]/gscomp[0]+(float)rand()/RMAX-0.5)][(short)(ge_time[i]/gscomp[2]+(float)rand()/RMAX-0.5)]++;
    }
  }
}



//------------------------------------------------------------------------------



//sort alpha-gated matrix
if(histogram_option==9){
  for(i=1;i<=mb_fold;i++){
    if(alphas>0){ 
      for(j=1;j<=clean_ge;j++){
        if(ge_energy[j]/gscomp[0] < MAT_LENGTH-1 && ge_energy[j]/gscomp[0] > 0){
          for(k=1;k<=clean_ge;k++){
            if(k!=j && ge_energy[k]/gscomp[0]>0 && ge_energy[k]/gscomp[0]<MAT_LENGTH-1){
              particle_number++;
              matrix[0][(short)(ge_energy[j]/gscomp[0]+(float)rand()/RMAX-0.5)][(short)(ge_energy[k]/gscomp[0]+(float)rand()/RMAX-0.5)]++;
            }
          }
        }
      }
    }
  }
}


//------------------------------------------------------------------------------



//sort microball-gated matrices
if(histogram_option==901){
  for(i=1;i<=mb_fold;i++){
    if(particle_counter[i]==1){ 
      for(j=1;j<=clean_ge;j++){
        if(ge_energy[j]/gscomp[0] < MAT_LENGTH-1 && ge_energy[j]/gscomp[0] > 0){
          for(k=1;k<=clean_ge;k++){
            if(k!=j && ge_energy[k]/gscomp[0]>0 && ge_energy[k]/gscomp[0]<MAT_LENGTH-1){
              particle_number++;
              matrix[0][(short)(ge_energy[j]/gscomp[0]+(float)rand()/RMAX-0.5)][(short)(ge_energy[k]/gscomp[0]+(float)rand()/RMAX-0.5)]++;
            }
          }
        }
      }
    }

    if(particle_counter[i]==2){
      for(j=1;j<=clean_ge;j++){
        if(ge_energy[j]/gscomp[0] < MAT_LENGTH-1 && ge_energy[j]/gscomp[0] > 0){
          for(k=1;k<=clean_ge;k++){
            if(k!=j && ge_energy[k]/gscomp[0]>0 && ge_energy[k]/gscomp[0]<MAT_LENGTH-1){
              particle_number++;
              matrix[1][(short)(ge_energy[j]/gscomp[0]+(float)rand()/RMAX-0.5)][(short)(ge_energy[k]/gscomp[0]+(float)rand()/RMAX-0.5)]++;
            }
          }
        }
      }
    }
  }
}



//------------------------------------------------------------------------------



// sort gated time vs ratio maps
if(histogram_option==902){ 

   for(i=1;i<=mb_fold;i++){
      if(mb_pass[i][2]/mbcomp[2] > 0 && mb_pass[i][2]/mbcomp[2] < MAT_LENGTH && mb_pass[i][5] > 0 &&
         mb_pass[i][5] < MAT_LENGTH && particle_counter[i]==1 || particle_counter[i]==2){
         matrix[mb_pass[i][4]][(short)(mb_pass[i][5]+(float)rand()/RMAX-0.5)][(short)(mb_pass[i][2]/mbcomp[2]+(float)rand()/RMAX-0.5)]++;
         particle_number++;
      }
   }
}



//------------------------------------------------------------------------------



// count particles
if(histogram_option==10){

  for(i=1;i<=mb_fold;i++){

    if(particle_counter[i]==1) spec[1][mb_pass[i][4]]++;
    if(particle_counter[i]==2) spec[2][mb_pass[i][4]]++;
      
    
  }
}



