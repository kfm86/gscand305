void writespec(char filename[40], int length, float psp[]);
void writemat(char filename[40], short mat_number);
short sidechannelcorrection(short ge_id_sc, short ge_energy, short ge_side, float offset, float gain, float thresh, float gscomp_sc);
float sctheta(short ge_id_sc, float ge_theta_sc);
int polygate(float (*polygon)[4], int *point);
void read_polygons(char polyfile[40], int npolygons, float (*polyarray)[20][4]);
void checkconflicts();
void read_diskfilelist(char diskfilelistfile[40], char diskfilename[100][60]);

