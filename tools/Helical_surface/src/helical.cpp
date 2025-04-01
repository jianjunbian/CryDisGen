/*==========================================================

  To create mesh grid of a single helical surface

	==========================================================*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

const double PI = 3.1415926;

typedef struct {
	double x, y, z;
} Vect;

double VDot(Vect a, Vect b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}
double VLen(Vect a)
{
	return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

int main(int argc, char **argv)
{
	Vect xg, yg, zg;     // global 
	Vect xl, yl, zl;     // local 

	if(argc != 28){
		printf("Usage:\n");
		printf("     Helical_surface -xg a b c -yg a b c -zg a b c ");
		printf("-xl a b c -yl a b c -zl a b c ");
		printf("r theta pitch\n");
		printf("\n");
		printf("Description:\n");
		printf("     To create the mesh grid of a helical surface\n");
		printf("\n");
		printf("Options:\n");
		printf("     -xg: global coordination system\n");
		printf("     -yg: \n");
		printf("     -zg: \n");
		printf("     -xl: local coordination system\n");
		printf("     -xl: \n");
		printf("     -xl: \n");
		printf("       r: radius of helical\n");
		printf("   theta: rotation angle (in degree)\n");
		printf("   pitch: pitch size of the helical surface\n");
		exit(1);
	}
  xg.x = atof(argv[2]);
  xg.y = atof(argv[3]);
  xg.z = atof(argv[4]);

  yg.x = atof(argv[6]);
  yg.y = atof(argv[7]);
  yg.z = atof(argv[8]);

  zg.x = atof(argv[10]);
  zg.y = atof(argv[11]);
  zg.z = atof(argv[12]);

  xl.x = atof(argv[14]);
  xl.y = atof(argv[15]);
  xl.z = atof(argv[16]);

  yl.x = atof(argv[18]);
  yl.y = atof(argv[19]);
  yl.z = atof(argv[20]);

  zl.x = atof(argv[22]);
  zl.y = atof(argv[23]);
  zl.z = atof(argv[24]);

	double r, theta, pitch;
	r  = atoi(argv[25]);
	theta  = atof(argv[26])/180 * PI;
	pitch  = atof(argv[27]);

	int nr = 8;    // mesh count along radial direction
	int nt = 40;   // mesh count along theta direction
  double dr = r / nr;
	double dtheta = theta / nt;
	double dz = pitch / 2. / PI;
	double zmax = theta * dz;
	printf("----------------------------------------\n");
	printf(   "Height of helical: %f\n", zmax);
	printf("----------------------------------------\n");
	
	double x, y, z;

	//----------------------------------------------------
	// Matrix for coordination transformation
	//----------------------------------------------------
	double m[3][3], np[3], op[3];
	m[0][0] = VDot(xl, xg)/VLen(xl)/VLen(xg);
	m[1][0] = VDot(xl, yg)/VLen(xl)/VLen(yg);
	m[2][0] = VDot(xl, zg)/VLen(xl)/VLen(zg);

	m[0][1] = VDot(yl, xg)/VLen(yl)/VLen(xg);
	m[1][1] = VDot(yl, yg)/VLen(yl)/VLen(yg);
	m[2][1] = VDot(yl, zg)/VLen(yl)/VLen(zg);

	m[0][2] = VDot(zl, xg)/VLen(zl)/VLen(xg);
	m[1][2] = VDot(zl, yg)/VLen(zl)/VLen(yg);
	m[2][2] = VDot(zl, zg)/VLen(zl)/VLen(zg);

	//----------------------------------------------------
	// array to store the mesh grid, which can be used to 
	// visualize the surface by ovito
	//----------------------------------------------------
  Vect pos[nr+1][nt+1];
	//----------------------------------------------------
	FILE *fp = fopen("helical.point.dump", "w");
	fprintf(fp, "ITEM: TIMESTEP\n");
	fprintf(fp, "0\n");
	fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
	fprintf(fp, "%d\n", (nr+1) * (nt+1));
	//fprintf(fp, "%d\n", natom - natom_del + ngridh[0]);
	fprintf(fp, "ITEM: BOX BOUNDS pp pp pp\n");
	fprintf(fp, "%f %f\n", -r, r);
	fprintf(fp, "%f %f\n", -r, r);
	fprintf(fp, "%f %f\n", -zmax * 0.5, zmax * 0.5);
	fprintf(fp, "ITEM: ATOMS id type x y z\n");

	for(int i = 0; i <= nr; i ++){
		for(int j = 0; j <= nt; j ++){
			x = dr * i * cos(dtheta * j);
			y = dr * i * sin(dtheta * j);
			z = dz * dtheta * j;
			z -= zmax / 2.;
			//-------------------------------------------------------------------------
			// coordination transform from local to global cartesian coordinate system
			//-------------------------------------------------------------------------
			for(int s = 0; s < 3; s ++) np[s] = 0;
			op[0] = x; op[1] = y; op[2] = z;

			for(int s = 0; s < 3; s ++){
				for(int t = 0; t < 3; t ++){
					np[s] += m[s][t] * op[t];
				}
			}
			x = np[0]; y = np[1]; z = np[2];
			//-------------------------------	 */
			pos[i][j].x = x;
			pos[i][j].y = y;
			pos[i][j].z = z;
			int id = 1;
			fprintf(fp, "%d %d %10f %10f %10f\n", id, 1, x, y, z);
		}
	}
	fclose(fp);

	//----------------------------------------------------
	// output the mesh grid used by CryDisGen
	//----------------------------------------------------
	fp = fopen("mesh.dat", "w");
	// construct triangulate mesh grid
	int v0, v1, v2, vn;
	int ntri = nr * nt * 2; // total triagulates
	fprintf(fp, "%d\n", ntri);

	for(int i = 0; i < nr; i ++){
		for(int j = 0; j < nt; j ++){
			fprintf(fp, "%10f %10f %10f\n", pos[i][j].x, pos[i][j].y, pos[i][j].z);
			fprintf(fp, "%10f %10f %10f\n", pos[i+1][j].x, pos[i+1][j].y, pos[i+1][j].z);
			fprintf(fp, "%10f %10f %10f\n", pos[i][j+1].x, pos[i][j+1].y, pos[i][j+1].z);

			fprintf(fp, "%10f %10f %10f\n", pos[i][j+1].x, pos[i][j+1].y, pos[i][j+1].z);
			fprintf(fp, "%10f %10f %10f\n", pos[i+1][j].x, pos[i+1][j].y, pos[i+1][j].z);
			fprintf(fp, "%10f %10f %10f\n", pos[i+1][j+1].x, pos[i+1][j+1].y, pos[i+1][j+1].z);
		}
	} 
	fclose(fp);

	return 0;
}
