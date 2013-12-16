#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/ioctl.h>
#include <sys/mtio.h>
#include <time.h>
#include <math.h>
#include "scandef.h"
#include "scanfunction.h"

FILE *fpr, *fp, *fps;

// MATRIXSIZE for short or int matrices toggle next two lines
//short matrix[NMAT][MAT_LENGTH][MAT_LENGTH];
int matrix[NMAT][MAT_LENGTH][MAT_LENGTH];
int     ndo[22], histogram_option;
int     disk_nblocks[300], disk_block_skip[300], ndiskfiles;
float   beta;
int main(int argc, char *argv[]){

   FILE     *fp_par, *fp_gs,*fp_ic,*fp_mb, *fp_stat, *fp_rng, *fp_mbc, *fp_v, *fp_dssd, *fp_dc, *fp_dl, *fp_sd, *fp_data;
   FILE     *fp_op, *fp_de, *fp_cal;
   int      fd, fd_out, i,j,k, mm,ll,block_skip, file_skip, nblocks, newblock, newblockflag;
   char     devname[60], output_devname[60], *ptr;
   char     ans1[1],ans2[1], ans3[1], ans4[1];
   char     makedir[100]="mkdir ";
   char     output_disk_file[100],temp_disk_file[100];
   short    getout;
   int      gs_det, gs_ring[111], mb_det, ge_id_sc, mb_fold;
   int      ic_lim1[50],ic_lim2[50],n_ic_cor;
   float    ic_de1_cor[50],ic_de2_cor[50],ic_de3_cor[50];
   float    ic_de_f,ic_de_total_f;
   int      ic_de_cor;
   int      counter[110];
   int      eff_number=0,old_eff_number, modeflags;
   float    mbcomp[3], gscomp[3];
   float    gs_theta[111], gs_phi[111], ge_theta_sc, ge_theta;
   float    mb_theta[96], mb_phi[96], ratio, beta;
   unsigned short    buff[8192];
   unsigned short    *p,*event_start, *event_end;
   char     gs_config_header[80], mb_config_header[80], config_start[1];
   char     parameter_file[40],gs_config_file[40],mb_config_file[40];
   char     icfile[40];
   char     scan_stat_file[40], vlines_file[40], dssdlookupfile[40], dssdconfigfile[40];
   char     scan_dir[80];
   char     test_dir[80], spec_file_temp[80], mat_file_temp[80];
   char     implantcalfile[40], dssdcalfile[40];
   struct   mtop mt_command;
   struct   mtget mt_status;
   time_t   start, finish, now;
   struct   tm *time_ptr;
   double   duration;
   float    doppler;

   unsigned short *op, obuff[8192],*start_of_output_event;
   int      oi,output_event_length;

   int      chan_p[11], chan_h[11];
   int      drop_protons, drop_alphas, max_files;
   int      last_word, last_ffff;
   short    block_length, header_length;
   unsigned short data_length, *ptr_data_length;
   int      event_length;
   int      count;
   int      clean_ge, dirty_ge, bgo_only, total_fold;
   int      usec_clock_l, usec_clock_m, usec_clock_h;
   double   usec_clock, usec_previous, usec_previous2, usec_start;
   double   usec_clock_blockstart, usec_clock_lastblockstart;
   double   usec_difference, usec_since, msec_since, sec_since;
   int      min_clock, int_clock, sec_clock, msec_clock;
   int      min_start, int_start, sec_start, msec_start;
   int      min_difference ,int_difference, sec_difference, msec_difference;
   int      usec_temp, usec_marker;
   int      tac1, tac2, rf, bgo_low_res_sum, ge_low_res_sum;
   int      bgo_hit_pattern[30], ge_hit_bit[30], ge_id[30],iid;
   int      ge_pileup[30], ge_energy[30], ge_side[30], ge_time[30];
   float    ge_e_f[30];
   int      ge_overflow[30], fold_overflow;
   int      total_fera_words;
   int      data0,data1, data2, data3, data4, data5, data6, data7, data8, data9, data10;
   int      data11, data12, data13, data14, data15;
   unsigned short *fera_start;
   int      vsn, fera_data, fera_words, fera_chan;
   int      mball_id, mball[102][6], mb_pass[102][6];
   int      mball_rank, mball_rank_fold[3];
   int      ppac_left, ppac_right, ppac_up, ppac_down, ppac_rf, ppac_dssd;
   int      rl,ud;
   int      ppac_x, ppac_y;
   int      ppac_de, grecoil_tac, ppac_rf_tof,rftof,rftof2, tgppac, ic_de1, ic_de2, ic_de3;
   int      tacppacrf, ppac_tof;
   int      ic_de1_uni, ic_de2_uni, ic_de3_uni;
   int      ic_de_total, ic_de_cath,ic_de, icwords, ic_de1_raw,ic_de2_raw,ic_de3_raw;
   int      ic_de_raw1,ic_de_raw2;
   int      sibe, sibtac,etsq1,etsq2,etsq3,etsq,tof,tof2;
   int      length;
   char     extension[4];
   float    inc_spec[SPEC_LENGTH];
   int      spec[NSPEC+1][SPEC_LENGTH];
   char     spec_file[80], mat_file[80];
   short    ii,jj,qq,nn;
   int      microball, fma, ionchamber;
   int      implantf, implantb, decayf, decayb, ppac, decay;
   short    ge_id_corrected;
   float    sc_gain[111], sc_offset[111], sc_thresh[111], time_offset1[111], time_offset2[111], time_factor[111], time_a;
   float    gain, offset, thresh;
   int      sccounter[250];
   int      increment;

   /* particle ID */
   int     protons, alphas, unidentified_particles, channel, number_particles, particle_test, particle_number;
   int     abs_type[97], q, absorber;
   float   abs_thickness[97], al_thickness[97];
   float   costh3[97], sinth3[97];
   float   thden, kpr, v3, Vcm, Ecm, Eocm, th3cm;
   float   beam_energy, tgt_thickness;
   int     tot_a, tgt_a, tgt_z, beam_a, beam_z;
   float   Eplab_MeV[10], Ealab_MeV[10], E_ap;
   int     p_csi[10], a_csi[10];
   float   Elab_MeV, particle_mass;
   float   sincm, sinlab, sigcm;
   int     csi_num;
   float   e3cm[10][10][10], phi3cm[10][10][10];
   int     Ep[10], csi_p[10], Ea[10], csi_a[10], particle_counter[20];

   /* gammasphere calibration coefficients */
   char    gs_cal_file[40];
   float   ge_coeff_a[120],ge_coeff_b[120],ge_coeff_c[120],ge_coeff_d[120],ge_coeff_e[120],ge_coeff_f[120];
   float   part1, part2, part3, part4, part5, part6;
   int     n_ge, id, order;
   char    gs_cal_header[80];

   /* neutrons */
   int     neutron[32][6], neutron_id, neutrons, neutron_fold;
   int     neutron_rank, neutron_rank_fold[6], ndet_pass[32][6];

   /* polygon stuff */
   char    polyfile[40];
   float   polygon[20][4];
   float   polyarray[POLYGONS+1][20][4];
   float   (*poly_ptr)[4];
   float   (*array_ptr)[20][4];
   int     point[2], npolygons, in=0, polynumber;
   int     n_vlines, v_index;
   float   vlo_x[97],vlo_y[97], vhi_x[97], vhi_y[97], vline[97][1024], vgrad;
   /* end of polygon stuff */

   /* microball calibration stuff */
   char    mball_cal_file[40];
   char  mball_cal_header[80];
   int   n_mb, mb_id;
   float mb_E_offset[97],mb_t_offset[97];
   float   slope[97], pedestal[97], acoeff[97], bcoeff[97], ccoeff[97];
   float   mbc1, mbc2, mbc3, mbc4, mbc5;

   /*
    * dssd
   unsigned int striprf[50], striprb[50], erf[50], erb[50];
   int   nstrips, dssdlookup[170] ,index;
   unsigned int stripdf[50], stripdb[50], edb[50], edf[50];
   int   silena_words, silena_pattern, silena_chan, silena_data;
   int   silena_word_count;
   int   ppac_sib_tof;
   int   dssdindex, implindex, dssdt;
   float    dssda[190], dssdb[190], dssdra[190], dssdrb[190], impla[100], implb[100];
   int     dssdflag[200],dssdthreshold[200];
   double   pixel[2][48][48];
   int      decayhitpattern, decayfold, decayenergy[200], decayid, decay_channel;
   int      decayfoldtotal, decayfoldtotalf, decayfoldtotalb;
   int      decayadc, n_export_adcs;
   int      recoilid, recoilenergy[200], recoilfold, recoil;
   int      recoilflag;
   int      silenaA[17],silenaB[17];
   int      recoileventid[200], recoileventenergy[200], recoileventxy[200], recoileventfold;
   int      decayeventid[200], decayeventenergy[200], decayeventxy[200], decayeventfold;
   int      decayenergydifference, decayeventenergymean;
   int      recoileventidfront[100], recoileventidback[100];
   int      recoileventenergyfront[100], recoileventenergyback[100];
   int      recoileventfoldfront, recoileventfoldback;
   int      decayeventidfront[100], decayeventidback[100];
   int      decayeventenergyfront[100], decayeventenergyback[100];
   int      decayeventfoldfront, decayeventfoldback;
   double   recoilarray[90][90][2];
   int      gamma_array[90][90][100];
   short    m;
   int      correlated_gamma_energy[100], correlated_gamma_fold;
   int      correlated_ppac_x;*/

   char     diskfilename[300][60], diskfilelistfile[40];
   int      diskfilecount, total_blocks;

   /* PREPROCESSING */

   if(argc<2){
      printf("Usage: %s parameter_file\n",argv[0]);
      exit(0);
   }

   strcpy(parameter_file,argv[1]);
   start=time(0);
   time(&now);
   time_ptr=localtime(&now);

   printf("\n---------------------------------------------");
   printf("\n >>> GSCAN305d - gs-mb sort");
   printf("\n >>> john f. smith - manchester - july 1998");
   printf("\n >>> version for gsfma305 - november 2012");
   printf("\n >>> also modified by kieran f. mulholland");
   printf("\n >>> comments: John.F.Smith@uws.ac.uk");
   printf("\n--------------------------------------------\n");

   printf("\nscan started: %s\n",asctime(time_ptr));

   if((fp_par = fopen(parameter_file, "r"))==NULL){
      printf("cannot open %s\n", parameter_file);
      exit(-1);
   }else  printf("opened file %s\n", parameter_file);

   printf("\n");

   /* initialize all counters */
   for(i=0;i<110;i++) counter[i]=0;
   for(i=0;i<250;i++) sccounter[i]=0;
   getout=0;

   /* --- initialise usec clock --- */
   usec_start=msec_start=sec_start=min_start=usec_marker=0;

   if(DISKSCAN){
      if(!fscanf(fp_par,"%s",diskfilelistfile)){
         printf("error reading disk filelist name: %s\n",diskfilelistfile);
         exit(-1);
      }else printf("Diskscan disk files: \t%s\n",diskfilelistfile);
   }
   if(PRESCAN){
      if(!fscanf(fp_par,"%s",output_devname)){
         printf("error reading prescan output device name: %s\n",output_devname);
         exit(-1);
      }else printf("Prescan output: \t\t%s\n",output_devname);
   }

   /*if(TAPESCAN){
      fscanf(fp_par,"%d",&block_skip);
      fscanf(fp_par,"%d",&file_skip);
      fscanf(fp_par,"%d",&nblocks);
      fscanf(fp_par,"%d",&max_files);
   }*/

   if(!fscanf(fp_par,"%s",gs_config_file)) exit(-1);
   if(!fscanf(fp_par,"%s",mb_config_file)) exit(-1);
   fscanf(fp_par,"%s",scan_stat_file);
   fscanf(fp_par,"%s",scan_dir);
   for(i=1; i<21; i++) fscanf(fp_par,"%d",&ndo[i]);
   fscanf(fp_par,"%d",&histogram_option);
   fscanf(fp_par,"%f %f %f", &gscomp[0], &gscomp[1], &gscomp[2]);
   fscanf(fp_par,"%f %f %f", &mbcomp[0], &mbcomp[1], &mbcomp[2]);
   if(ndo[1]) printf("mbcomp[1],[2],[3] = %f %f %f\n",mbcomp[0],mbcomp[1],mbcomp[2]);
   fscanf(fp_par,"%f",&beta);
   fscanf(fp_par,"%s",gs_cal_file);
   fscanf(fp_par,"%s",mball_cal_file);
   fscanf(fp_par,"%s",polyfile);

   printf("the following ndo logicals are set:\t");
   for(i=1; i<21; i++)
      if(ndo[i])
         printf("\t%d",i);
   printf("\n");
   printf("histogram option %d has been chosen\n", histogram_option);

   if(TAPESCAN){
      if(nblocks == -1){
         printf("\nno block limit specified\n");
         nblocks = 300000;
      } else printf("\nscanning %d blocks\n",nblocks);
   }

   /* --- Diskscan initialise --- */
   if(DISKSCAN) read_diskfilelist(diskfilelistfile,diskfilename);
   checkconflicts();
   if(DISKSCAN) max_files=ndiskfiles;

   if(max_files < 0) max_files=1000;
   if(max_files <500){
      printf("scanning %d file",max_files);
      if(max_files == 1) printf("\n");
      else printf("s\n");
   }
   else printf("scanning all files\n");

   /* --- Prepare scan directory --- */
   strcpy(test_dir,scan_dir);
   strcat(test_dir,"/scan_info.dat");
   if((fp_sd = fopen(test_dir, "w"))==NULL){
           printf("\n\nThe scan directory\n%s\ndoes not exist...",scan_dir);
           strcat(makedir,scan_dir);
           system(makedir);
      strcpy(makedir,"mkdir ");
   }

   if((fp_sd = fopen(test_dir, "w"))==NULL){
           printf("\n\nError error error: - cannot open scan directory:\n");
           printf("%s\n", scan_dir);
           printf("(Please check that it exists)\n\n");
           exit(-1);
   }
   fclose(fp_sd);

   /* --- Open config/stat files --- */
   if((fp_gs = fopen(gs_config_file, "r"))==NULL){
      printf("cannot open %s\n", gs_config_file);
      exit(-1);
   }else printf("opened file %s\n", gs_config_file);

   if((fp_mb = fopen(mb_config_file, "r"))==NULL) {
      printf("cannot open %s\n", mb_config_file);
      exit(-1);
   }else printf("opened file %s\n", mb_config_file);

   if((fp_stat = fopen(scan_stat_file, "w"))==NULL){
      printf("cannot open %s\n", scan_stat_file);
      exit(-1);
   }else printf("opened file %s\n", scan_stat_file);

   /* --- Print info to screen --- */
   fprintf(fp_stat,"\nSCAN INFORMATION\n");
   time(&now); time_ptr=localtime(&now);
   fprintf(fp_stat,"Scan started at: %s",asctime(time_ptr));
   fprintf(fp_stat,"Histogram option: %d\n",histogram_option);
   fprintf(fp_stat,"Disk file list: %s\n",diskfilelistfile);
   fprintf(fp_stat,"Scan directory: %s\n",scan_dir);
   fprintf(fp_stat,"The following ndo flags are set: ");
   for(i=1; i<21; i++){
      if(ndo[i]) fprintf(fp_stat,"%d  ",i);
   }
   fprintf(fp_stat,"\n\n");

   if(TAPESCAN){
      if((fd = open(devname, O_RDONLY))==-1){
         printf("cannot open tape drive %s\n",devname);
         exit(1);
      }else printf("tape opened for reading\n");
   }

   /* --- Open the output tape for the prescan --- */
   if(PRESCAN && OUTPUT_TAPE){
      if((fd_out = open(output_devname, O_RDWR))==-1){
         printf("cannot open output tape drive\n");
         exit(1);
      }
   }else printf("output tape opened for writing\n");
   printf("\n");

   /* --- create the output directory for the prescan --- */
   if(PRESCAN && OUTPUT_DISK){
      strcpy(test_dir,output_devname);
      strcat(test_dir,"/prescan_info.dat");
      if((fp_sd = fopen(test_dir, "w"))==NULL){
         printf("\n\nThe prescan output directory\n");
         printf("%s\n", output_devname);
         printf("does not exist...");
         strcat(makedir,output_devname);
         system(makedir);
         strcpy(makedir,"mkdir ");
         printf(":created\n");
      }

      sprintf(extension,"%d",counter[82]);
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
   } // End prescan output disk setup

   /* --- GS config setup --- */
   for(i=0;i<20;i++)
      if(!fgets(gs_config_header,80,fp_gs)) printf("error with gammasphere config file header");

   for(i=1; i<111; i++){
      fscanf(fp_gs,"%d %d %f %f %f %f %f %f %f",&gs_det, &gs_ring[i], &gs_theta[i], &gs_phi[i],
            &sc_gain[i], &sc_offset[i], &sc_thresh[i], &time_offset1[i], &time_offset2[i]);
      time_offset1[i]*=2;
      time_offset2[i]*=2;
      time_factor[i]=(4740-4000)/(time_offset2[i]-time_offset1[i]);
      if(ndo[1])
         printf("\t%d \t%d \t%f \t%f \t%f \t%f \t%f \t%f \t%f\n",gs_det, gs_ring[i], gs_theta[i],
               gs_phi[i], sc_gain[i], sc_offset[i], sc_thresh[i], time_offset1[i], time_offset2[i]);
   }

   /* --- uB config setup --- */
   if(!fgets(mb_config_header,70,fp_mb)) printf("error with mb file header");
   printf("reading the microball config file\n\n");
   for(i=1; i<97; i++){
      fscanf(fp_mb,"%d %f %f %f %d %f",&mb_det, &mb_theta[i], &mb_phi[i], &al_thickness[i],
            &abs_type[i], &abs_thickness[i]);
      costh3[i]=cos(mb_theta[i]*3.1415927/180.0);
      sinth3[i]=sin(mb_theta[i]*3.1415927/180.0);
      printf("%d %f %f %f %d %f %f %f\n",mb_det, mb_theta[i], mb_phi[i], al_thickness[i],
            abs_type[i], abs_thickness[i], sinth3[i], costh3[i]);
   }

   fclose(fp_mb);
   fclose(fp_par);
   fclose(fp_gs);

   // Read the Gammasphere calibration coefficients file
   if((fp_cal = fopen(gs_cal_file, "r"))==NULL){
      printf("cannot open %s\n", gs_cal_file);
      exit(-1);
   }else printf("opened file %s\n", gs_cal_file);

   n_ge=0;

   for(i=0;i<10;i++){
      if(!fgets(gs_cal_header,80,fp_cal)) printf("error with gammasphere config file header");
      printf("%s",gs_cal_header);
   }

   fscanf(fp_cal,"%d",&n_ge);
   printf("Number of germanium detectors: %d\n",n_ge);

   for(i=1; i<=n_ge; i++){
      fscanf(fp_cal,"%d %d %e %e %e %e %e %e",&id,&order,&ge_coeff_a[i],&ge_coeff_b[i],
            &ge_coeff_c[i],&ge_coeff_d[i],&ge_coeff_e[i],&ge_coeff_f[i]);
      ge_coeff_a[i]=(float)ge_coeff_a[i]*3.0;
      ge_coeff_b[i]=(float)ge_coeff_b[i];
      ge_coeff_c[i]=(float)ge_coeff_c[i]/3;
      ge_coeff_d[i]=(float)ge_coeff_d[i]/9;
      ge_coeff_e[i]=(float)ge_coeff_e[i]/27;
      ge_coeff_f[i]=(float)ge_coeff_f[i]/81;
      printf("%d %d \n%e %e %e \n%e %e %e\n",id,order,ge_coeff_a[i],ge_coeff_b[i],ge_coeff_c[i],ge_coeff_d[i],ge_coeff_e[i],ge_coeff_f[i]);
   }
   fclose(fp_cal);
   // End of: Read the Gammasphere calibration coefficients file

   // Read the Microball calibration coefficients file
   if((fp_mbc = fopen(mball_cal_file, "r"))==NULL){
      printf("cannot open %s\n", mball_cal_file);
      exit(-1);
   }else printf("opened file %s\n", mball_cal_file);

   for(i=0;i<10;i++){
      if(!fgets(mball_cal_header,80,fp_mbc))printf("error with microball config file header");
      printf("%s",mball_cal_header);
   }
   n_mb=0;
   fscanf(fp_mbc,"%d",&n_mb);
   printf("Number of Microball detectors: %d\n",n_mb);

   for(i=1;i<n_mb;i++){
      fscanf(fp_mbc,"%d %e %e %e %e %e %f\n",&mb_id,&mb_E_offset[i],&mb_t_offset[i],&acoeff[i],&bcoeff[i],&ccoeff[i],&pedestal[i]);
      mb_E_offset[i]=(float)mb_E_offset[i]*2;
      mb_t_offset[i]=(float)mb_t_offset[i]*2;
      acoeff[i]=acoeff[i]/2;
      bcoeff[i]=bcoeff[i];
      ccoeff[i]=ccoeff[i]*2;
      pedestal[i]=pedestal[i]*2;
      printf("%d %e %e %e %e %e %f\n",mb_id,mb_E_offset[i],mb_t_offset[i],acoeff[i],bcoeff[i],ccoeff[i],pedestal[i]);
   }
   fclose(fp_mbc);
   // End of: Read the Microball calibration coefficients file

   array_ptr=polyarray;
   read_polygons(polyfile, POLYGONS, array_ptr);

   diskfilecount=-1;
   disksortstart: ;
   if(DISKSCAN){
      diskfilecount++;
      counter[3]++;
      nblocks=disk_nblocks[diskfilecount];
      counter[4]=counter[5]=counter[6]=0;
      if((fp_data = fopen(diskfilename[diskfilecount], "rb"))==NULL){
         printf("ERROR! cannot open %s\n", diskfilename[diskfilecount]);
         goto disksortstart;
         exit(-1);
      }
      else printf("\nopened file %s\n", diskfilename[diskfilecount]);
   }

   /* PROCESSING */

   start: ;
   if(ans1[0] == 'x'){
      while(ans2[0] != 'r'){
         printf("Put new input tape in the drive and enter r when ready\n");
         scanf("%s",ans2);
         printf("Waiting 40 seconds to make sure tape is ready\n");
         sleep(40); /* wait 40 seconds to make sure that tape is ready */
         if((fd = open(devname, O_RDONLY))==-1){
       printf("cannot open tape drive\n");
       exit(1);
         }else printf("tape opened for reading\n");
         file_skip=0; block_skip=0;
      }
      strncpy(ans2,"x",sizeof(ans2));
   }

   /* --- determine number of blocks to sort --- */
   if(nblocks == -1){
      printf("\nno block limit specified\n");
      nblocks = 1000000;
   }else printf("\nscanning %d blocks\n",nblocks);

   /* file skip must come before block skip */
   if(TAPESCAN && file_skip){
      printf("skipping %d file",file_skip);
      if(file_skip != 1) printf("s\n");
      else printf("\n");
      mt_command.mt_op = MTFSF;
      mt_command.mt_count = file_skip;
      ioctl(fd, MTIOCTOP, &mt_command);
   }
   if(TAPESCAN && block_skip){
      printf("skipping %d blocks\n",block_skip);
      mt_command.mt_op = MTFSR;
      mt_command.mt_count = block_skip;
      ioctl(fd, MTIOCTOP, &mt_command);
   }
   /* --- Read events into buffer --- */
   while(counter[5] < 6 && counter[3] <= max_files && counter[6] <= nblocks){
      if(TAPESCAN) i=read(fd, buff, MAX_BYTES);
      if(DISKSCAN) i=fread(buff,1,16384,fp_data);
      counter[1]++;

      if(i != 16384){
         if(i == 112){
            printf("block size = 112\n");
            counter[2]++;
            continue;
         }
         else if(i == 90){
            printf("block size = 90\n");
            counter[3]++;
            counter[0]++;
            printf("FILE HEADER BLOCK\n");
            continue;
         }
         else if(i == 0){
            printf("block size = 0\n");
            counter[5]++;
            continue;
         }
         else{
            printf("block size = something else");
            counter[4]++;
            printf("INCORRECT BLOCK LENGTH (%d)",i);
            continue;
         }
      } // End block size check

      /* if you haven't read a file-header within the first 20 blocks then you must have started mid file,
         therefore increment file-header counter anyway */
      if(counter[6] == 20 && counter[3] == 0) counter[3] = 1;

      /* now have a data buffer - extract events from it */
      if(SWAPBYTES){
         for(i=0;i<8192;i++){
            if(ndo[17]) printf("BEFORE: \tbuff[%d]=%x\n",i,buff[i]);
            buff[i]=((buff[i]&0x00ff)<<8)|((buff[i]&0xff00)>>8);
            if(ndo[17]) printf("AFTER: \tbuff[%d] %x\n",i,buff[i]);
         }
      }

      last_ffff=8191;

      while(buff[last_ffff] != 0xffff) last_ffff--;

      if(last_ffff <= 1) continue; /* no ffff in the block */

      counter[6]++;
      newblockflag=1;

      if(fmod(counter[6],500) == 0){
         time(&now);
         time_ptr=localtime(&now);
         printf("file: %d - block number: %d -- %s", counter[3],counter[6],asctime(time_ptr));
      }

      if(fmod(counter[6],10000) == 0){
         time(&now);
         time_ptr=localtime(&now);
         fprintf(fp_stat,"file: %d - block number: %d -- %s", counter[3],counter[6],asctime(time_ptr));
      }

      p = buff;

   /* read useful words from the block header */

      old_eff_number=eff_number;
      block_length=*(p+1);
      eff_number=*(p+4);
      if(ndo[1]){printf("eff number: %d blocks: %d\n", eff_number, counter[6]);}
      header_length=*(p+3)/2;
      modeflags=*(p+7);
      data_length=*(p+8);
      ptr_data_length = &data_length;

      newblock=1;
      p+=header_length;

      if(ndo[1] == 1){
         printf("data blocks = %d",counter[6]);
         printf("header_length = %d, last_ffff= %d\n", header_length,last_ffff);
         printf("p = %hn, data_length = %d, ptr_data_length = %d\n", p, data_length, *ptr_data_length);
      }

      if(header_length != 11) continue;

      while(p-buff < data_length){
         if(ndo[1]){printf("Start of next event p-buff %ld\n",p-buff);}
         if((*p & 0xf000) != 0x8000 || *(p+(*p & 0x00ff)-1) != 0xffff){
            counter[16]++;
            if(ndo[1]){printf("Counter[16] incremented HERE \n");}
            while(*p != 0xffff && p-buff<data_length){p++;}
            p++;
            continue;
         }

         if(ndo[1]){printf("Middle of good event checks\n");}

         if(*p == 0xffff){
            p++;
            counter[7]++;
            continue;
         }

         if(ndo[1]){printf("Before good event checks\n");}

         /* have a good event */
         counter[8]++;
         microball=0;
         fma=0;
         ionchamber=0;
         implantf=0;
         implantb=0;
         decayf=0;
         decayb=0;
         ppac=0;
         etsq=tac1=tac2=rf=0;
         event_length = *p & 0x00ff;
         event_start = p;
         event_end =  p+event_length-1;

         if(ndo[1]){printf("event_length = %d, event_start = %hn,event_end = %d\n\n",
               event_length, event_start, *(event_end));}

         p++;
         clean_ge = *p & 0x00ff;
         if(clean_ge > 30){
            fold_overflow++;
            while(*p != 0xffff) p++;
            continue;
         }

         p++;
         dirty_ge = (*p & 0xff00)>>8;
         bgo_only = *p & 0x00ff;
         total_fold=clean_ge+dirty_ge+bgo_only;

         p++;
         usec_clock_h = *p;

         p++;
         usec_clock_m = *p;

         p++;
         usec_clock_l = *p;
         usec_previous2=usec_previous;
         usec_previous=usec_clock;

         if(newblockflag) usec_clock_lastblockstart=usec_clock_blockstart;
         usec_clock=pow(2,31)*usec_clock_h+pow(2,16)*usec_clock_m+usec_clock_l;

         if(ndo[1]){
            printf("usec_clock/500 s = %f\n",usec_clock/500e6);
            printf("usec_clock_h = %d\nusec_clock_m = %d\nusec_clock_l = %d\n"
                  ,usec_clock_h,usec_clock_m,usec_clock_l);
            printf("usec_clock = %e\n",usec_clock);
            printf("counter[8] = %d\n",counter[8]);
         }

         if(counter[8]==1){
            if(usec_start==0){usec_start=usec_clock;}
         }

      /* --- usec clock restarts after run 16. marker to avoid overlapping usec clock values--- */
         if(usec_clock<usec_start&&usec_marker==0){usec_marker++;}
         if(usec_marker==1){usec_clock+=4e11;} //note 4e11 is arbritary value

         if(newblockflag){usec_clock_blockstart=usec_clock;}

         msec_clock=(int)(usec_clock/1e3);
         sec_clock=(int)(usec_clock/1e6);
         min_clock=(int)(usec_clock/60e6);
         int_clock=min_clock&0x7ff;

         if(counter[8]==1){
            if(msec_start==0) msec_start=msec_clock;
            if(sec_start==0) sec_start=sec_clock;
            if(min_start==0) min_start=min_clock;
         }

         usec_difference=usec_clock-usec_start;
         msec_difference=msec_clock-msec_start;
         sec_difference=sec_clock-sec_start;
         min_difference=min_clock-min_start;

         if(ndo[1]){
            printf("\nusec_difference=%e %e %e\n",usec_difference, usec_start, usec_clock);
            printf("int %d\n",(int)usec_difference);
            printf("msec_clock %d msec_start %d msec_difference %d\n",msec_clock, msec_start, msec_difference);
            printf("sec_clock %d sec_start %d sec_difference %d\n",sec_clock, sec_start, sec_difference);
            printf("min_clock %d min_start %d min_difference %d\n",min_clock, min_start, min_difference);
            if(counter[8]>10000) exit(1);
            printf("usec clock value: %f diff: %f\n", usec_clock, usec_clock-usec_previous);
         }

         p++;
         tac1 = *p & 0x0fff;

         p++;
         rf = *p & 0x0fff;

         p++;
         ge_low_res_sum = *p;

         p++;
         bgo_low_res_sum = *p;

         if(ndo[1]){
            printf("\n----------NEW EVENT--------------\n");
            printf("clean_ge \t\t%d\n",clean_ge);
            printf("dirty_ge \t\t%d\n",dirty_ge);
            printf("bgo_only \t\t%d\n",bgo_only);
            printf("tac1 \t\t%d\n",tac1);
            printf("rf \t\t%d\n",rf);
            printf("microsec_clock_h \t\t%d\n",usec_clock_h);
            printf("microsec_clock_m \t\t%d\n",usec_clock_m);
            printf("microsec_clock_l \t\t%d\n",usec_clock_l);
            printf("\t min_clock%d int_clock %d\n", min_clock,int_clock); /* line added 13/07/04 */
            printf("ge_low_res_sum \t\t%d\n",ge_low_res_sum);
            printf("bgo_low_res_sum \t\t%d\n",bgo_low_res_sum);
         }

         for(i=1; i<=100;i++){
            for(j=1;j<5;j++){
               mball[i][j]=0;
            }
         }
         for(i=1;i<31;i++){
            for(j=1;j<6;j++){
               neutron[i][j]=0;
            }
         }

         neutrons=0;
         // --- end of event header --- //

         // --- start of data words --- //
         for(i=1; i<=clean_ge; i++){
            p++;
            bgo_hit_pattern[i] = (*p & 0xfe00)>>9;
            ge_hit_bit[i] = (*p & 0x0100)>>8;
            ge_id[i] = *p & 0x00ff;

            if(ge_id[i] > 110 || ge_id[i] < 0){
               ge_id[i]=111;
            }

            p++;
            ge_pileup[i] = (*p & 0x8000)>>15;
            ge_overflow[i] = (*p & 0x4000)>>14;
            ge_energy[i] = *p & 0x3fff;
            p++;
            ge_side[i] = *p & 0x0fff;
            p++;
            ge_time[i] = (*p & 0x1fff);

            if(ndo[1]){
               printf("ge_id[i] %d time_offset %d\n",ge_id[i],ge_time[i]);
               printf("---- %d ----\n", i);
               printf("bgo hit pattern %d\n", bgo_hit_pattern[i]);
               printf("ge hit bit %d\n", ge_hit_bit[i]);
               printf("ge id %d\n", ge_id[i]);
               printf("ge pileup %d\n", ge_pileup[i]);
               printf("ge overflow %d\n", ge_overflow[i]);
               printf("ge energy %d\n", ge_energy[i]);
               printf("ge side %d\n", ge_side[i]);
               printf("ge time %d\n", ge_time[i]);
            }
         }      // --- end of data words --- //

         // --- Start: germanium checks and repacking --- //

         if(ndo[2]){ // --- Align germanium times --- //
            for(i=1;i<=clean_ge;i++){
               //ge_time[i]=ge_time[i]+(4000-time_offset1[ge_id[i]]);
               time_a=ge_time[i]-time_offset1[ge_id[i]];
               time_a=(time_a*time_factor[ge_id[i]])+(float)rand()/RAND_MAX-0.5;
               ge_time[i]=time_a+4000;
            }
         } // --- End ndo[2] --- //

         if(ndo[7]){ // --- Gain match the Gammasphere energies for detectors that need it--- //
            for(i=1; i<=clean_ge; i++){
               if(fmod(ge_id[i],1000)==61||fmod(ge_id[i],1000)==71||fmod(ge_id[i],1000)==90||
                     fmod(ge_id[i],1000)==92||fmod(ge_id[i],1000)==106){
                  //printf("Before: ge_id: %d ge_energy: %d\n",ge_id[i],ge_energy[i]);
                  part1=ge_coeff_a[ge_id[i]];
                  //printf("part1 %f\n",part1);
                  part2=ge_coeff_b[ge_id[i]]*ge_energy[i]+(float)rand()/RAND_MAX-0.5;
                  //printf("part2 %f\n",part2);
                  part3=ge_coeff_c[ge_id[i]]*ge_energy[i]*ge_energy[i]+(float)rand()/RAND_MAX-0.5;
                  //printf("part3 %f\n",part3);
                  part4=ge_coeff_d[ge_id[i]]*ge_energy[i]*ge_energy[i]*ge_energy[i]+(float)rand()/RAND_MAX-0.5;
                  //printf("part4 %f\n",part4);
                  part5=ge_coeff_e[ge_id[i]]*ge_energy[i]*ge_energy[i]*ge_energy[i]*ge_energy[i]+(float)rand()/RAND_MAX-0.5;
                  //printf("part5 %f\n",part5);
                  part6=ge_coeff_f[ge_id[i]]*ge_energy[i]*ge_energy[i]*ge_energy[i]*ge_energy[i]*ge_energy[i]+(float)rand()/RAND_MAX-0.5;
                  //printf("part6 %f\n",part6);
                  ge_energy[i]=part1+part2+part3+part4+part5+part6;
                  if(ge_energy[i]<0){
                     ge_energy[i]=0;
                  }
                  //printf("After: ge_id: %d ge_energy: %d\n",ge_id[i],ge_energy[i]);
               }
            }
         } // --- End ndo[7]

         if(ndo[8]){ // --- Before time-gating checks --- //
            if(ndo[1]){
               printf("\n BEFORE REPACKING EVENT NUMBER: %d\n", counter[8]);
               printf("\nClean germanium fold: %d\n", clean_ge);
               for(i=1;i<=clean_ge;i++){
                  printf("\nGERMANIUM NUMBER: %d\n", i);
                  printf("ge energy: \t%d\n", ge_energy[i]);
                  printf("ge time: \t%d\n", ge_time[i]);
                  printf("ge side: \t%d\n", ge_side[i]);
                  printf("ge id: \t%d\n", ge_id[i]);
                  printf("ge hit bit: \t%d\n", ge_hit_bit[i]);
                  printf("ge pileup: \t%d\n", ge_pileup[i]);
                  printf("ge overflow: \t%d\n", ge_overflow[i]);
                  printf("bgo hit pattern: \t%d\n", bgo_hit_pattern[i]);
               }
            }
         } // -- End before time-gating checks

         if(ndo[8]==1){ // --- basic time gating --- //
            for(i=1;i<=clean_ge;i++){
               if(ge_time[i]<2*1998 || ge_time[i]>2*2017){
                  for(j=i;j<clean_ge;j++){
                     ge_id[j]=ge_id[j+1];
                     ge_energy[j]=ge_energy[j+1];
                     ge_time[j]=ge_time[j+1];
                     ge_side[j]=ge_side[j+1];
                     ge_overflow[j]=ge_overflow[j+1];
                     ge_pileup[j]=ge_pileup[j+1];
                     ge_hit_bit[j]=ge_hit_bit[j+1];
                     bgo_hit_pattern[j]=bgo_hit_pattern[j+1];
                  }
                  clean_ge--;
                  i--;
               }
            }
            total_fold=clean_ge+dirty_ge+bgo_only; // Redefine total fold in case of time gating
         } // --- end ndo[8]==1 --- //

         if(ndo[8]==2){ // --- polygonal time gating --- //
            for(i=1;i<=clean_ge;i++){
               point[0]=(int)ge_energy[i]/gscomp[0];
               poly_ptr=polyarray[1];
               point[1]=(int)ge_time[i]/gscomp[2];

               if(ndo[1]){
                  printf("\npoint[0]: \t%d\n", point[0]);
                  printf("point[1]: \t%d\n", point[1]);
               }

               if(!polygate(poly_ptr, point)){
                  if(ndo[1]){
                     printf("bad ge_energy[i]: %d\n", ge_energy[i]);
                     printf("bad ge_time[i]: %d\n", ge_time[i]);
                  }
                  for(j=i;j<clean_ge;j++){
                     ge_id[j]=ge_id[j+1];
                     ge_energy[j]=ge_energy[j+1];
                     ge_time[j]=ge_time[j+1];
                     ge_side[j]=ge_side[j+1];
                     ge_overflow[j]=ge_overflow[j+1];
                     ge_pileup[j]=ge_pileup[j+1];
                     ge_hit_bit[j]=ge_hit_bit[j+1];
                     bgo_hit_pattern[j]=bgo_hit_pattern[j+1];
                  }
                  clean_ge--;
                  i--;
               }
            }
            total_fold=clean_ge+dirty_ge+bgo_only; // Redefine total fold in case of time gating
         } // --- end ndo[8]==2 --- //

         if(ndo[8]){ // After time-gating checks
            if(ndo[1]){
               printf("\n AFTER REPACKING EVENT NUMBER: %d\n", counter[8]);
               printf("\nClean germanium fold: %d\n", clean_ge);
               for(i=1;i<=clean_ge;i++){
                  printf("\nGERMANIUM NUMBER: %d\n", i);
                  printf("ge energy: \t%d\n", ge_energy[i]);
                  printf("ge time: \t%d\n", ge_time[i]);
                  printf("ge side: \t%d\n", ge_side[i]);
                  printf("ge id: \t%d\n", ge_id[i]);
                  printf("ge hit bit: \t%d\n", ge_hit_bit[i]);
                  printf("ge pileup: \t%d\n", ge_pileup[i]);
                  printf("ge overflow: \t%d\n", ge_overflow[i]);
                  printf("bgo hit pattern: \t%d\n", bgo_hit_pattern[i]);
               }
            }
         } // --- End after time-gating checks --- //

         if(ndo[3]){    /* --- add side channel information to the detector numbers --- */
            for(i=1; i<=clean_ge; i++){
               if(sc_thresh[ge_id[i]] != 4000.0 && ge_energy[i] > 100){
                  sccounter[ge_id[i]]++;
                  if(ge_side[i]==0){
                     sccounter[ge_id[i]+120]++;
                  }
                  gain=sc_gain[ge_id[i]];
                  offset=sc_offset[ge_id[i]];
                  thresh=sc_thresh[ge_id[i]];
                  ge_id_corrected=sidechannelcorrection(ge_id[i], ge_energy[i], ge_side[i], offset, gain, thresh, gscomp[1]);
                  ge_id[i]=ge_id_corrected;
               }
               else ge_id[i]+=1000;
            }
         } // --- end ndo[3] --- //

         if(ndo[20]){
            if(fmod(counter[6],1000)==0){
               for(i=1;i<111;i++){
                  printf("\t%d %d %d\n",i,sccounter[i],sccounter[120+i]);
               }
            }
         } // --- end ndo[20] --- //

         if(ndo[1]){
            for(i=1;i<=clean_ge;i++){
               printf("BEFORE i %d ge_id[i] %d ge_energy[i] %d\n",i, ge_id[i], ge_energy[i]);
            }
         }

         if(ndo[4]){ // Doppler correction - requires ndo[3]
            for(i=1; i<=clean_ge;i++){
               if(ge_id[i]-fmod(ge_id[i],1000) != 1000 && ge_id[i]-fmod(ge_id[i],1000) != 2000){
                  ge_id_sc = ge_id[i];
                  ge_theta_sc = gs_theta[(int)fmod(ge_id[i],1000)];
                  ge_theta = sctheta(ge_id_sc, ge_theta_sc);
               }else{
                  ge_theta = gs_theta[(int)fmod(ge_id[i],1000)];
               }

               /* then Doppler correct germanium energies*/
               if(ndo[1]){printf("\t%d %d %f %f\n",i, ge_energy[i], ge_theta_sc, ge_theta);}

               doppler=(1-beta*cos(ge_theta*3.1415927/180))/sqrt(1-pow(beta,2));
               if(ndo[1]){printf("doppler=%f\n",doppler);}

               /*ge_energy[i]=(int)( (float)(ge_energy[i] + (float)rand()/32768.0-0.5) * doppler);*/
               ge_energy[i]=(int)((float)(ge_energy[i] + (float)rand()/RAND_MAX-0.5) * doppler);

               if(ndo[1]){printf("after Doppler correction %d\n",ge_energy[i]);}
            }
         } // --- end ndo[4] --- //

         if(ndo[5]){ // Doppler correct WITHOUT side channel info - requires != ndo[3]
            for(i=1; i<=clean_ge;i++){
               ge_theta = gs_theta[(int)fmod(ge_id[i],1000)];
               doppler=(1-beta*cos(ge_theta*3.1415927/180));
               doppler=(1-beta*cos(ge_theta*3.1415927/180))/sqrt(1-pow(beta,2));

               if(ndo[1]){printf("doppler=%f\n",doppler);}

               //ge_energy[i]=(int)( (float)(ge_energy[i] + (float)rand()/32768.0-0.5) * doppler);
               ge_energy[i]=(int)( (float)(ge_energy[i] + (float)rand()/RAND_MAX-0.5) * doppler);

               if(ndo[1]){printf("after Doppler correction %d\n",ge_energy[i]);}
            }
         } // --- end ndo[5] --- //

         if(ndo[1]){
            for(i=1;i<=clean_ge;i++){
               printf("AFTER i %d ge_id[i] %d ge_energy[i] %d\n",i, ge_id[i], ge_energy[i]);
            }
         }
         p++;

         if(!GSONLY){
            if(*p != 0xffff){
               if(*p != 0xff00){
                  if(ndo[1]){
                     printf("\nmissing ff00\n");
                     printf("\n\ndata blocks = %d\n",counter[6]);
                  }
                  continue;
               }

               p++;
               total_fera_words = *p;

               if(ndo[1]){printf("total_fera_words %d\n",total_fera_words);}

               fera_start = p;

               while(p < fera_start + total_fera_words){
                  p++;
                  if((*p & 0x8000) == 0x8000){
                     //this is a fera header word
                     vsn = *p & 0x00ff;
                     if(vsn < 0x51 || vsn > 0x5c){
                        fera_words = (*p & 0x7800)>>11;
                        if(fera_words == 0){fera_words = 16;}
                        if(ndo[1]){
                           printf("vsn \t\t\t%d %x\n", vsn, vsn);
                           printf("fera_words \t\t\t%d\n", fera_words);
                           printf("data blocks = %d", counter[6]);
                        }

                        // now read the fera data words

                        for(j=1; j<= fera_words; j++){
                           p++;
                           if((*p & 0x8000)>>15 != 0){continue;}
                           fera_chan = (*p & 0x7800)>>11;
                           fera_data = *p & 0x07ff;
                           if(ndo[1]){
                              printf("fera chan  = %d\n", fera_chan);
                              printf("fera data  = %d\n", fera_data);
                           }

//--------------------------------------------------------------------------

                           if((vsn & 0x00f0)==0x0060||(vsn & 0x00f0)==0x0070||(vsn & 0x00f0)==0x0080){
                              if((vsn & 0x000f) >= 1 && (vsn & 0x000f) <= 6){
                                 mball_id = ((vsn & 0x000f)-1)*16 + fera_chan + 1;
                                    if(ndo[1] == 1){printf("mball_id %d\n",mball_id);}
                                    if(!microball){microball++;}
                                    if((vsn & 0x00f0) == 0x0060){ // this is a microball energy
                                       mball[mball_id][1] = fera_data;
                                    }
                                    if((vsn & 0x00f0) == 0x0070){ // this is a microball time
                                       mball[mball_id][2] = fera_data-rf+2200;
                                    }
                                    if((vsn & 0x00f0) == 0x0080){ // this is a microball pid
                                       mball[mball_id][3] = fera_data;
                                    }
                              }
                           } // end of reading mb data
                        } // end of for(j<=fera_words) loop
                     } // end of vsn check
                  } // end of fera header word (*p&0x8000==0x8000) check
               }  // end of fera words loop


//--------------------------------------------------------------------------

               if(microball && !IGNORE_MICROBALL){
                  mball_rank_fold[0]=0;
                  mball_rank_fold[1]=0;
                  mball_rank_fold[2]=0;
                  mb_fold=0;
                  for(i=0;i<97;i++){
                     for(j=1;j<=4;j++){
                        mb_pass[i][j]=0;
                     }
                  }
                  for(i=0;i<20;i++){particle_counter[i]=0;}

                  for(i=0;i<97;i++){
                     mball_rank=0;
                     if(mball[i][1] != 0 && mball[i][1] > pedestal[i]){mball_rank++;}
                     if(mball[i][2] != 0){mball_rank++;}
                     if(mball[i][3] != 0){mball_rank++;}
                     if(mball_rank==3){
                        mball_rank_fold[2]++;
                        mb_pass[mball_rank_fold[2]][1] = mball[i][1]; /* energy */
                        mb_pass[mball_rank_fold[2]][2] = mball[i][2]; /* time */
                        mb_pass[mball_rank_fold[2]][3] = mball[i][3]; /* pid */
                        mb_pass[mball_rank_fold[2]][4] = i;           /* detector id */

                        //some detectors have much lower gain
                        if(mb_pass[mball_rank_fold[2]][4]==20){
                           mb_pass[mball_rank_fold[2]][1]=
                                 2*(mb_pass[mball_rank_fold[2]][1]+((float)rand()/RAND_MAX)-0.5);
                           mb_pass[mball_rank_fold[2]][3]=
                                 2*(mb_pass[mball_rank_fold[2]][3]+((float)rand()/RAND_MAX)-0.5);
                        }

                        //some detectors have much higher gain
                        if(mb_pass[mball_rank_fold[2]][4]==38 || mb_pass[mball_rank_fold[2]][4]==31 ||
                              mb_pass[mball_rank_fold[2]][4]==34 || mb_pass[mball_rank_fold[2]][4]==48
                              || mb_pass[mball_rank_fold[2]][4]==59){

                           mb_pass[mball_rank_fold[2]][1]=
                                 (mb_pass[mball_rank_fold[2]][1]+((float)rand()/RAND_MAX)-0.5)/2;
                           mb_pass[mball_rank_fold[2]][3]=
                                 (mb_pass[mball_rank_fold[2]][3]+((float)rand()/RAND_MAX)-0.5)/2;
                        }

                        if(ndo[6]){ // --- Calibrate Microball --- //
                           mb_pass[mball_rank_fold[2]][1]+=mb_E_offset[i]; /* energy */
                           mb_pass[mball_rank_fold[2]][2]+=mb_t_offset[i]; /* time */
                           mb_pass[mball_rank_fold[2]][3]-=(acoeff[i]*mb_pass[mball_rank_fold[2]][1]*mb_pass[mball_rank_fold[2]][1]); /* pid */
                           if(ccoeff[i]>=0){
                              mb_pass[mball_rank_fold[2]][3]-=ccoeff[i]; /* move pid down*/
                           }else{
                              mb_pass[mball_rank_fold[2]][1]+=ccoeff[i]; /* move pid left */
                           }
                           if(mb_pass[mball_rank_fold[2]][4]==26){
                              mb_pass[mball_rank_fold[2]][3]-=28;
                           }
                           //line up uB times after gain shift
                           if(usec_clock>7.565e11) mb_pass[mball_rank_fold[2]][2]+=96;

                           //cut out noise in time signal
                           if(mb_pass[mball_rank_fold[2]][4]==1 || mb_pass[mball_rank_fold[2]][4]==6 ||
                              mb_pass[mball_rank_fold[2]][4]==7 || mb_pass[mball_rank_fold[2]][4]==27 ||
                              mb_pass[mball_rank_fold[2]][4]==66 || mb_pass[mball_rank_fold[2]][4]==69){
                              if(mb_pass[mball_rank_fold[2]][2]<1360 || mb_pass[mball_rank_fold[2]][2]>2200){
                                 if(usec_clock>5.222976e10 && usec_clock < 7.979232e10){
                                    mb_pass[mball_rank_fold[2]][1]=0;
                                    mb_pass[mball_rank_fold[2]][2]=0;
                                    mb_pass[mball_rank_fold[2]][3]=0;
                                    mball_rank--;
                                    mball_rank_fold[2]--;
                                 }
                              }
                           }
                        } // --- End ndo[6] --- //

                        if(ndo[1]){
                           printf("mb ID: %d\nmb E: %d\nmb t: %d\nmb PID: %d\n",
                              mb_pass[mball_rank_fold[2]][4], mb_pass[mball_rank_fold[2]][1],
                              mb_pass[mball_rank_fold[2]][2], mb_pass[mball_rank_fold[2]][3]);
                        }

                        ratio=120*((mb_pass[mball_rank_fold[2]][1]+((float)rand()/RAND_MAX)-0.5)/(mb_pass[mball_rank_fold[2]][3]));
                        mb_pass[mball_rank_fold[2]][5] = (int)ratio;
                        mb_fold=mball_rank_fold[2];
                     }//end of microball_rank == 3 check

                     if(mball_rank == 2) mball_rank_fold[1]++;
                     if(mball_rank == 1) mball_rank_fold[0]++;
                  } //end of microball data-read loop

                  if(ndo[9]){  /* --- particle gating --- */
                     protons=0;
                     alphas=0;
                     unidentified_particles=0;
                     number_particles=0;

                     for(i=1;i<=mb_fold;i++){
                        particle_test=0;
                        point[0]=mb_pass[i][5];
                        point[1]=(int)(mb_pass[i][2]/mbcomp[2]);

                        poly_ptr=polyarray[mb_pass[i][4]+1]; /* --- alpha gate --- */
                        if(polygate(poly_ptr,point)){
                           alphas++;
                           number_particles++;
                           particle_counter[i]=1;
                           Ea[alphas]=mb_pass[i][1];
                           csi_a[alphas]=mb_pass[i][4];
                        }else if(!(polygate(poly_ptr,point))){
                           particle_test+=100;
                        }else{
                           printf("This shouldn't happen - bogus particle ID");
                           printf("mb ID: %d\nmb E: %d\nmb t: %d\nmb PID: %d\n",
                                 mb_pass[mball_rank_fold[2]][4], mb_pass[mball_rank_fold[2]][1],
                                 mb_pass[mball_rank_fold[2]][2], mb_pass[mball_rank_fold[2]][3]);
                        }

                        poly_ptr=polyarray[mb_pass[i][4]+100]; /* --- proton gate --- */
                        if(polygate(poly_ptr,point)){
                           protons++;
                           number_particles++;
                           particle_counter[i]=2;
                           Ep[protons]=mb_pass[i][1];
                           csi_p[protons]=mb_pass[i][4];
                        }else if(!(polygate(poly_ptr,point))){
                           particle_test+=100;
                        }else{
                           printf("This shouldn't happen - bogus particle ID");
                           printf("mb ID: %d\nmb E: %d\nmb t: %d\nmb PID: %d\n",
                                 mb_pass[mball_rank_fold[2]][4], mb_pass[mball_rank_fold[2]][1],
                                 mb_pass[mball_rank_fold[2]][2], mb_pass[mball_rank_fold[2]][3]);
                        }

                        if(ndo[1]){
                           printf("particle gating:\n");
                           printf(" alphas = %d\n protons = %d\n number_particles = %d\n",
                                 alphas,protons,number_particles);
                           printf(" particle_test = %d\n mb_fold = %d\n i = %d\n",
                                 particle_test,mb_fold,i);
                        }
                     }

                     if(protons==0 && alphas==0) unidentified_particles++;

                     counter[20]+=protons;
                     counter[21]+=alphas;
                     counter[22]+=unidentified_particles;
                     counter[30]+=number_particles;
                  } /* --- end ndo[9] --- */

                  if(ndo[1]){
                     printf("\t\t0: %d\n",mball_rank_fold[0]);
                     printf("\t\t1: %d\n",mball_rank_fold[1]);
                     printf("\t\t2: %d\n",mball_rank_fold[2]);
                     for(i=1;i<=mball_rank_fold[2];i++) printf(":::: %d\n",mb_pass[i][4]);
                     for(k=1;k<=96;k++){
                        printf("%d\t%d\t%d\t%d\n",k,mball[k][1],mball[k][2],mball[k][3]);
                     }
                     printf("\n\ndata_blocks = %d\n",counter[6]);
                  }

               } /* end of microball */

               newblockflag=0;

            } // end of if(p!=0xffff)
         } // end of if fera data
         else p=event_end;

         p++;

         #include "histogram_options_gsfma305.c"
         //#include "prescan_option_gsfma305.c"

      }
/* --------------------------- end of event ------------------------------ */

   }
/* --------------------------- end of files ------------------------------ */

   total_blocks+=counter[6];
   printf("\n----------------------------------------------------\n");
   printf("End of file %d. Blocks sorted in this file %d\n",counter[3],counter[6]);
   printf("Total blocks sorted: %d\n",total_blocks);
   printf("----------------------------------------------------\n\n");
   fprintf(fp_stat,"\n----------------------------------------------------\n");
   fprintf(fp_stat,"End of file %d. Blocks sorted in this file %d\n",counter[3],counter[6]);
   fprintf(fp_stat,"Total blocks sorted: %d\n",total_blocks);
   fprintf(fp_stat,"----------------------------------------------------\n\n");

   fclose(fp_data);

   if(DISKSCAN && counter[3] < max_files){
      goto disksortstart;
   }

   /* --- POSTPROCESSING --- */

   counter[31]=particle_number; /* --- total number of particle-gated matrix increments in hist option 9** --- */

   if(counter[5]>=5 || counter[3]>=max_files || counter[6]>=nblocks && !getout){
      /* counters: [5]=consecutive zero length blocks [3]=files sorted [6]=blocks sorted */
      printf("\nEnd of tape or specified amount of data has been read\n\n");

      if(PRESCAN){
         ans1[0]='j';
         while(ans1[0] != 'y' && ans1[0] != 'Y' && ans1[0] != 'n' && ans1[0] != 'N'){
            printf("Do you want to continue with another input file or tape?\n");
            scanf("%s",ans1);
         }
         if(ans1[0] == 'y' || ans1[0] == 'Y'){
            getout=0; close(fd);
            strcpy(ans1,"x");
            counter[50]+=counter[3];
            counter[51]+=counter[5];
            counter[52]+=counter[6];
            counter[3]=counter[4]=counter[5]=counter[6]=0;
         }
         else if(ans1[0] == 'n' || ans1[0] == 'N'){
            getout=1; close(fd);
            mt_command.mt_op = MTWEOF;
            mt_command.mt_count=1;
            ioctl(fd_out, MTIOCTOP, &mt_command);
            close(fd_out);
         }
      }
      else getout=1;
   }

   if(!getout) goto start;
   if(getout){
      printf("\n\ndata_blocks = %d\n",counter[6]);
      printf("\n\nfile count = %d\n",counter[0]);
      printf("\n\nzero_length_blocks = %d\n",counter[5]);
      printf("\n\ndata_length = %d\n",data_length);
      printf("\n\nmax files = %d\n",max_files);

      if(SPECTRA){
         for(i=0;i<=NSPEC;i++){
            for(j=1;j<=SPEC_LENGTH;j++){
               inc_spec[j]=0;
            }
            strcpy(spec_file,SPEC_FILE);
            printf("\tspectrum # %d\n",i);
            sprintf(extension,"%d",i);
            strcat(spec_file,extension);
            strcat(spec_file,".spe");
            strcpy(spec_file_temp,spec_file);
            strncpy(spec_file,scan_dir,sizeof(spec_file));
            strcat(spec_file,"/");
            strcat(spec_file,spec_file_temp);
            for(j=1;j<=SPEC_LENGTH;j++){
               inc_spec[j]=(float)(spec[i][j]);
            }
            writespec(spec_file,SPEC_LENGTH,inc_spec);
         }
      }

      fps=fopen("spec1temp.dat", "w");
      for(i=1;i<2048;i++){
         fprintf(fps,"%d %d\n",i,spec[1][i]);
      }
      fclose(fps);

      printf("test002\n");

      if(MATRICES){
         for(i=0;i<NMAT;i++){
            strcpy(mat_file,MAT_FILE);
            printf("\tmatrix # %d\n",i);
            sprintf(extension,"%d",i);
            strcat(mat_file,extension);
            // MATRIXSIZE for short or int matrices toggle next two lines
            strcat(mat_file,".m4b");
            //strcat(mat_file,".mat");
            strcpy(mat_file_temp,mat_file);
            strcpy(mat_file,scan_dir);
            strcat(mat_file,"/");
            strcat(mat_file,mat_file_temp);
            writemat(mat_file,i);
         }
      }

      printf("test003\n");
      finish=time(0);
      duration=difftime(finish,start);
      printf("test004\n");
      fprintf(fp_stat,"\n\tScan statistics\n\n");
      for(i=1; i<101; i++){
         printf("\tcounter[%d] = %d\n",i,counter[i]);
         fprintf(fp_stat,"counter[%d] = %d\n",i,counter[i]);
      }

      printf("\n");
      time(&now);
      time_ptr=localtime(&now);
      printf("\nScan finished successfully at: %s\n",asctime(time_ptr));
      fprintf(fp_stat,"\nScan finished successfully at: %s\n",asctime(time_ptr));
      printf("Using histogram option %d\n\n",histogram_option);
      printf("(Scan duration: %f seconds)\n\n",duration);
      fprintf(fp_stat,"(Scan duration: %f seconds)\n\n",duration);
      fclose(fp_stat);
   }
} // end of main


void writespec(char filename[40], int length, float psp[]){

   FILE *fpr;
   char filename_cut[9];
   int ii=1, jj=24, kk, i;

   fpr=fopen(filename,"wb");
   kk=4*length;
   strncpy(filename_cut,filename,8);

   fwrite(&jj,4,1,fpr);
   fwrite(filename_cut,8,1,fpr);
   fwrite(&length,4,1,fpr);
   fwrite(&ii,4,1,fpr);
   fwrite(&ii,4,1,fpr);
   fwrite(&ii,4,1,fpr);
   fwrite(&jj,4,1,fpr);
   fwrite(&kk,4,1,fpr);
   fwrite(psp,4*length,1,fpr);
   fwrite(&kk,4,1,fpr);

   fclose(fpr);
}

void writemat(char filename[40], short mat_number){

   FILE *fp;
   int i,j;
   short *p;

   fp=fopen(filename,"wb");

   if(!fp){
      fprintf(stderr, "can't create file %s\n",filename);
      exit(1);
   }
   // MATRIXSIZE for short or int matrices toggle next two lines
   fwrite(matrix[mat_number],sizeof(int),MAT_LENGTH*MAT_LENGTH,fp);
   //fwrite(matrix[mat_number],sizeof(short),MAT_LENGTH*MAT_LENGTH,fp);

   fclose(fp);
}

short sidechannelcorrection(short ge_id_sc, short ge_energy, short ge_side, float sc_offset, float sc_gain, float sc_thresh, float gscomp_sc){
   short side_index, ge_id_corrected;
   short ge_energy_keV, ge_side_keV;
   short flag;

   flag=0;

   if(ndo[1]){
      printf("Side Channel Correction for Id=%d\n",ge_id_sc);
   }
   ge_energy_keV=(ge_energy)/3;
   if(ge_energy_keV < 1){
      ge_energy_keV=1;
   }
   ge_side_keV=ge_side*sc_gain+sc_offset;
   if(ge_side_keV < 1){
      ge_side_keV=1;
   }
   side_index=100;

   if(ge_energy_keV < 2*sc_thresh){
      side_index=10;
      if(ndo[1]){
         printf("total energy < 2*sc_thresh\n");
      }
   }
   else if(ge_energy_keV < 10*sc_thresh && ge_energy_keV >= 2*sc_thresh){
      if((float)ge_side_keV/ge_energy_keV < 0.5){
    side_index=21;
    if(ndo[1]){
            printf("side_index=21\n");
         }
      }
      else{
    side_index=22;
    if(ndo[1]){
            printf("side_index=22\n");
         }
      }
   }
   else if(ge_energy_keV >= 10*sc_thresh){
      if(ndo[1]){
         printf("(float) side/total = %f\n",(float)ge_side_keV/ge_energy_keV);
      }
      if((float)ge_side_keV/ge_energy_keV < 0.1){
    side_index=31;
    if(ndo[1]){
            printf("side_index=31\n");
         }
      }
      else if((float)ge_side_keV/ge_energy_keV > 0.9){
    side_index=32;
    if(ndo[1]){
            printf("side_index=32\n");
         }
      }
      else{
    side_index=10;
    if(ndo[1]){
            printf("side_index=10\n");
         }
      }
   }
   else{
      side_index=10;
      if(ndo[1]){
         printf("reached the end side index = 10\n");
      }
   }

   if(ndo[1]){
      printf("%d %d %d %d %d\n",ge_side_keV, ge_energy_keV,ge_side, ge_energy, (int)(1000*((float)ge_side_keV/ge_energy_keV)));
   }


   /* if(ge_id_sc==25||ge_id_sc==27||ge_id_sc==33||ge_id_sc==35||ge_id_sc==71||ge_id_sc==76)*/
   if(ge_id_sc==27 || ge_id_sc==34 || ge_id_sc==36 || ge_id_sc==38 || ge_id_sc==40 || ge_id_sc==71 || ge_id_sc==85){
      flag=1;
      if(ndo[1]){
         printf("Reflange detector! Id=%d\n,Flag=%d\n", ge_id_sc, flag);
      }
   }

   switch(side_index){
      case 100:
    printf("side_index = 100 SHOULD NOT HAPPEN");
    ge_id_corrected=ge_id_sc+1000;
    break;
      case 10:
         ge_id_corrected=ge_id_sc+2000;
         break;
      case 21:
         if(flag==0){
            ge_id_corrected=ge_id_sc+3000;
         }
         if(flag==1){
            ge_id_corrected=ge_id_sc+4000;
         }
         break;
      case 22:
         if(flag==0){
            ge_id_corrected=ge_id_sc+4000;
         }
         if(flag==1){
            ge_id_corrected=ge_id_sc+3000;
         }
         break;
      case 31:
         if(flag==0){
            ge_id_corrected=ge_id_sc+5000;
         }
         if(flag==1){
            ge_id_corrected=ge_id_sc+6000;
         }
         break;
      case 32:
         if(flag==0){
            ge_id_corrected=ge_id_sc+6000;
         }
         if(flag==1){
            ge_id_corrected=ge_id_sc+5000;
         }
         break;
   }
   return(ge_id_corrected);
} // --- end side channel correction --- //

float sctheta(short ge_id_sc, float ge_theta_sc){
   unsigned short location,flag,id;
   location=ge_id_sc-fmod(ge_id_sc,1000);

   switch(location){
      case 1000:
         printf("bogus - should not happen\n");
         break;
      case 2000:
         printf("bogus - should not happen\n");
         break;
      case 3000:
         ge_theta_sc=ge_theta_sc-2.4;
         break;
      case 4000:
         ge_theta_sc=ge_theta_sc+2.4;
         break;
      case 5000:
         ge_theta_sc=ge_theta_sc-2.9;
         break;
      case 6000:
         ge_theta_sc=ge_theta_sc+2.9;
         break;
   }
   return(ge_theta_sc);
}

void read_polygons(char polyfile[40], int npolygons, float (*polyarray)[20][4]){
  FILE *fp;
  int i=0, j=0, k;

  if ((fp = fopen(polyfile, "r"))==NULL) {
    printf("cannot open %s\n", polyfile);
    exit(-1);
  } else
    printf("opened file %s\n", polyfile);

  fscanf(fp,"%d",&k);
  npolygons=k;

  for(j=1;j<=npolygons;j++)
    {
      polyarray[j][19][0]=90;
      polyarray[j][19][1]=0;
      polyarray[j][19][2]=100000;
      polyarray[j][19][3]=0;

      fscanf(fp,"%f",&polyarray[j][18][0]); /* polygon index */
      fscanf(fp,"%f",&polyarray[j][18][1]); /* number of vertices */
      if(ndo[9])printf("read polygon %d\n",(int)polyarray[j][18][0]);

      for(i=0; i<(int)polyarray[j][18][1]; i++)
   {
     fscanf(fp,"%f %f",&polyarray[j][i][0],&polyarray[j][i][1]);
     if(ndo[9])printf("%d (%f %f)\n",(int)polyarray[j][18][0], polyarray[j][i][0], polyarray[j][i][1]);
     polyarray[j][i][2]=180*atan((double)polyarray[j][i][1]/polyarray[j][i][0])/3.1415927;
     polyarray[j][i][3]=sqrt(pow(polyarray[j][i][0],2)+pow(polyarray[j][i][1],2));
     if(polyarray[j][i][2] < polyarray[j][19][0]) polyarray[j][19][0]=polyarray[j][i][2];
     if(polyarray[j][i][2] > polyarray[j][19][1]) polyarray[j][19][1]=polyarray[j][i][2];
     if(polyarray[j][i][3] < polyarray[j][19][2]) polyarray[j][19][2]=polyarray[j][i][3];
     if(polyarray[j][i][3] > polyarray[j][19][3]) polyarray[j][19][3]=polyarray[j][i][3];
   }
      /* close the polygon (last vertex = first vertex) */
      for(i=0;i<5;i++)
   polyarray[j][(int)polyarray[j][18][1]][i]=polyarray[j][0][i];

    }
}

int polygate(float (*polygon)[4], int *point){

/* Explanation of polygate function. This function reads in data
   about a polygonal gate and a data point (x,y) and checks whether
   the data point is inside the polygon. If the point is inside the
   polygon the function returns 1 else it returns 0. The integer
   vertices is the number of vertices in the gate. The (x,y) coordinates
   of the data point are passed as (point[0],point[1]). The polygon array
   is as follows:
   polygon[i][0] = x coordinate of vertex i of polygon
   polygon[i][1] = y coordinate of vertex i of polygon
   polygon[i][2] = angle between position vector of vertex i
                   and the positive x axis counter-clockwise
   polygon[i][3] = distance from vertex i to the origin
   polygon[19][0] = minimum theta of any polygon vertex
   polygon[19][1] = maximum theta of any polygon vertex
   polygon[19][2] = minimum distance of any polygon vertex
   polygon[19][3] = maximum distance of any polygon vertex
   polygon[18][0] = polygon index
   polygon[18][1] = number of vertices MUST BE LESS THAN 18
   jfs january 1999*/

   float dist, angle;
   float xcp, ycp, distcp;
   float c, m1, m2;
   int cross[20], cross_count=0, in=0, i;

   for(i=0;i<20;i++){
      cross[i]=0;
   }

   dist=sqrt(pow(point[0],2)+pow(point[1],2));
   angle=180*atan((double)point[1]/point[0])/3.1415927;

   if(angle >= polygon[19][0] && angle <= polygon[19][1] && dist >= polygon[19][2] && dist <= polygon[19][3]){
      for(i=0;i<(int)polygon[18][1];i++){
    if(polygon[i][2] >= angle){
       if(polygon[i+1][2] < angle){
               cross[i]=1;
            }
         }
    if(polygon[i][2] <= angle){
            if(polygon[i+1][2] > angle){
               cross[i]=1;
            }
         }
      }

      for(i=0;i<(int)polygon[18][1];i++){
    if(cross[i]){
         /* next two lines to avoid infinite gradients */
       if(polygon[i][0]==polygon[i+1][0]){
               polygon[i][0]=1.01*polygon[i+1][0];
            }
       c=(polygon[i+1][0]*polygon[i][1]-polygon[i][0]*polygon[i+1][1])/(polygon[i+1][0]-polygon[i][0]);
       m1= (float) point[1]/point[0];
       m2= (polygon[i][1]-polygon[i+1][1])/(polygon[i][0]-polygon[i+1][0]);

       xcp=c/(m1-m2);
       ycp=m1*xcp;
       distcp=sqrt(pow(xcp,2)+pow(ycp,2));
       if(dist >= distcp){
               cross_count++;
            }
    }
      }

      if(cross_count != 0 && fmod(cross_count,2) != 0){
         return 1;
      }
      else{
         return 0;
      }
   }
   else{
      return 0;
   }
}


void checkconflicts(){

   if(DISKSCAN && !TAPESCAN){
      printf("Scanning data from disk...\n");
   }
   else if(!DISKSCAN && TAPESCAN){
      printf("Scanning data from tape...\n");
   }
   else if(TAPESCAN && DISKSCAN){
      printf("ERROR! Both disk and tape scan defined. Exiting.\n");
      exit(1);
   }
   else if(!TAPESCAN && !DISKSCAN){
      printf("ERROR! Neither disk nor tape scan defined. Exiting.\n");
      exit(1);
   }
   /*
   if(DISKSCAN && file_skip>0){
      printf("ERROR! Scanning data from disk and file-skip is non-zero..\n");
      exit(1);
   } */
}


void read_diskfilelist(char diskfilelistfile[40], char diskfilename[300][60]){

   FILE *fp;
   int i;

   if((fp = fopen(diskfilelistfile, "r"))==NULL){
      printf("cannot open %s\n", diskfilelistfile);
      exit(-1);
   }
   else{
      printf("opened file %s\n", diskfilelistfile);
   }

   fscanf(fp,"%d",&ndiskfiles);

   if(ndiskfiles>199){
      printf("Too many data files (maximum is 200).\n");
      exit(1);
   }
   printf("\nThe number of files to be sorted is %d\n\n",ndiskfiles);
   printf("The following files will be sorted...\n");
   for(i=0;i<ndiskfiles;i++){
      fscanf(fp,"%s %d %d",diskfilename[i], &disk_nblocks[i], &disk_block_skip[i]);
      printf("%d: \t%s (%d) (%d) \n",i, diskfilename[i], disk_nblocks[i], disk_block_skip[i]);
   }
   printf("\n");
}

/*
counters
0        files
1        total blocks
2        tape_header_blocks
3        file_header_blocks
4        bad_length_blocks
5        zero_length_blocks
6        16k data blocks
7        events that do not end 0xffff
8        good events
9        events with microball words
10       events with fma words
11       events with ion-chamber words
12       events with microball AND fma
13       events with microball AND fma AND ion-chamber
14       events with microball AND NOT fma
15       events with fma AND NOT microball
16       events that do not begin 0x8000 or end 0xffff
*/

