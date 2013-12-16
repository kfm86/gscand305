if(PRESCAN){
  
  /* create output buffer if this is first pass */
  if(counter[80]==0){
    for(i=0;i<OUT_BYTES/2;i++) obuff[i]=0;
    op=obuff;
    *op=0xbbbb; /* block header word */
    counter[80]++;
  }

  if(ndo[1]) printf("prescan_option position 1\n");
  op++;

  start_of_output_event=op;
  *op=0xe0e0;  /* will be overwritten later by event length */
  op++;

  if(ndo[1]) printf("prescan_option position 2\n");
  *op=(dirty_ge+bgo_only+clean_ge)<<8;
  *op=*op|clean_ge;
  op++;

  *op=ge_low_res_sum;
  op++;

  if(ndo[1]) printf("prescan_option position 3\n");
  *op=bgo_low_res_sum;
  op++;

  *op=rf;
  op++;

  if(microball) 
    *op=*op|(1<<15);

  if(ndo[1]) printf("prescan_option position 4\n");
  if(microball && protons < 15) 
    *op=*op|(protons<<4);

  if(ndo[1]) printf("prescan_option position 5\n");
  if(microball && alphas < 15) 
    *op=*op|alphas;

  if(ndo[1]) printf("prescan_option position 6\n");
  for(i=1;i<=clean_ge;i++){
    op++;
    *op=fmod(ge_id[i],1000);
    op++;
    *op=ge_energy[i];
  }
  if(ndo[1]) printf("prescan_option position 7.0\n");
  counter[109]=0;
  if(microball){
    for(i=1;i<=protons;i++){
      if(ndo[1]) printf("prescan_option position 7.1\n");
      op++;

      if(ndo[1]) printf("prescan_option position 7.2\n");
      *op=csi_p[i];
      op++;

      if(ndo[1]) printf("prescan_option position 7.3\n");
      *op=Ep[i];
      if(ndo[1]){
        printf("prescan_option position 7.4\n");
        counter[109]++; 
        if(counter[109]>10) exit(-1);
      }
    }
    for(i=1;i<=alphas;i++){
      if(ndo[1]) printf("prescan_option position 7.5\n");
      op++;

      if(ndo[1]) printf("prescan_option position 7.6\n");
      *op=csi_a[i];
      op++;

      if(ndo[1]) printf("prescan_option position 7.7\n");
      *op=Ea[i];
      
    }
  }
  if(ndo[1]) printf("prescan_option position 8\n");
  op++;

  *op=0xffff; /* --- end of event --- */

  output_event_length=op-start_of_output_event+1; /* work out event length including 0xffff */
  *(op-output_event_length+1)=output_event_length; /* --- make the first word of the event equal to the event length --- */

  if(op-obuff >= 8192-100-output_event_length){

    if(OUTPUT_DISK){
      oi=fwrite(obuff,1,OUT_BYTES,fp_op);
    }
    if(oi!=OUT_BYTES){
      printf("Problem writing output: is the output medium full?\n");
      getout=1;
      continue;
    }

    for(i=0;i<OUT_BYTES/2;i++) obuff[i]=0;
    
    op=obuff;
    *op=0xbbbb; /* block header word */
    counter[81]++;

    if(ndo[1]) printf("prescan_option position 10\n");

    if(counter[81]>0 && fmod(counter[81],OUTPUT_FILE_SIZE)==0){
    
      fclose(fp_op);
      counter[82]++;

      sprintf(extension,"%d",(int)counter[82]);
      strcpy(output_disk_file,DISK_FILE_STEM);
      strcat(output_disk_file,extension);
      strcat(output_disk_file,".dat");
      strcpy(temp_disk_file,output_disk_file);
      strcpy(output_disk_file,output_devname);
      strcat(output_disk_file,"/");  
      strcat(output_disk_file,temp_disk_file);  
    
      if((fp_op = fopen(output_disk_file, "wb"))==NULL){
        printf("cannot open %s\n", output_disk_file);
        exit(-1);
      }else printf("opened file %s\n", output_disk_file);

      counter[84]+=counter[81]; 
      counter[81]=0;
    
    }

    if(counter[81]>0 && fmod(counter[81],200)==0) printf("\t-> %d blocks written\n", counter[81]);
    
  }

  if(ndo[1]) printf("prescan_option position 11\n");

  if(counter[82] >= OUTPUT_FILE_NUMBER || counter[81]>=290000){
    printf("%d blocks written on this tape\n",counter[81]);
    mt_command.mt_op = MTWEOF;
    mt_command.mt_count=1;
    ioctl(fd, MTIOCTOP, &mt_command);
    printf("eom written and tape closed\n");
    close(fd_out);

    if(ndo[1]) printf("prescan_option position 12\n");

    while(ans4[0] != 'y' && ans4[0] != 'Y' && ans4[0] != 'n' && ans4[0] != 'N'){
      printf("Do you want to continue with another output tape?\n");
      scanf("%s",ans4);
    }

    if(ndo[1]) printf("prescan_option position 12.1\n");

    if(ans4[0] == 'y' || ans4[0] == 'Y'){
      ans3[0]='j';
      while(ans3[0] != 'r'){
        printf("Place a new output tape in the drive and enter r when ready\n");
        scanf("%s",ans3);
      }
          
      if((fd_out = open(output_devname, O_RDWR))==-1){
        printf("cannot open output tape drive\n");
        exit(1);
      }else printf("output tape opened for writing\n");
          
      if(ndo[1]) printf("prescan_option position 13\n");
          
      counter[81]=0;
      counter[82]=0;
      counter[83]++;
      strcpy(ans4,"x");
    }
    else if(ans4[0] == 'n' || ans4[0] == 'N'){
      getout=1; 
      break;
    }
  }			  
  if(ndo[1]) printf("prescan_option position 14\n");

} /* end of PRESCAN */


