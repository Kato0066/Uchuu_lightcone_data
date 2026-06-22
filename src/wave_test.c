#include <stdio.h>
#include <math.h>

#define N (128)

int main(int argc, char **arg){
  int i, j, k;
  double nu_cfl = 0.5;
  double c = 1.0;
  double dx = 1.0 / (double)N;
  double dt = nu_cfl*dx/c;
  double PI = 3.141592653589793;

  double tend;
  int istep=0;

  static double U0[N][N][N];
  static double U1[N][N][N];
  static double U2[N][N][N];

  //初期条件
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      for(k=0;k<N;k++)
        double x = (i+0.5)*dx;
      double y = (j+0.5)*dx;
      double z = (k+0.5)*dx;
      U0[i][j]= exp(-100.0*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));
      //U0[i][j]= sin(2.0*PI*x);
    }
  }


  double phi[N][N][N], A[N][N][N];
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){
      for(k=0; k<N; k++){
	double x = (i+0.5)*dx;
	double y = (j+0.5)*dx;
	double z = (k+0.5)*dx;
        phi[i][j][k] = 0.0;
        A[i][j][k] = 1.0/(1.0-4.0*phi[i][j][k])*(nu_cfl*nu_cfl);

	int ip = i+1 > N-1 ? i+1-N : i+1;
	int im = i-1 < 0 ? i-1+N : i-1;
	int jp = j+1 > N-1 ? j+1-N : j+1;
	int jm = j-1 < 0 ? j-1+N : j-1;
	int kp = k+1 > N-1 ? k+1-N : k+1;
	int km = k-1<0 ? k-1+N : k-1;
	U1[i][j][k] = U0[i][j][k] + 0.5 * A[i][j][k] * (U0[ip][j][k] + U0[im][j][k] + U0[i][jp][k] + U0[i][jm][k] + U0[i][j][kp] + U0[i][j][km] - 6.0 * U0[i][j][k]);
	//U1[i] = U0[i];
      }
    }

    double peak_positions[N];
    double times[N];
    double velocities[N - 1];

    double tnow = 0.0;
    while(tnow < tend) {
      for(i=0;i<N;i++){
	for(j=0;j<N;j++){
	  int ip = i+1 > N-1 ? i+1-N : i+1;
	  int im = i-1 < 0 ? i-1+N : i-1;
	  int jp = j+1 > N-1 ? j+1-N : j+1;
	  int jm = j-1 < 0 ? j-1+N : j-1;
	  int kp = k+1 > N-1 ? k+1-N : k+1;
	  int km = k-1<0 ? k-1+N : k-1;
	  U2[i][j][k] = 2.0 * U1[i][j][k] - U0[i][j][k] + A[i][j][k] * (U0[ip][j][k] + U0[im][j][k] + U0[i][jp][k] + U0[i][jm][k] + U0[i][j][kp] + U0[i][j][km] - 6.0 * U0[i][j][k]);
	}
      }

      for(i=0;i<N;i++){
	for(j=0;j<N;j++){
	  U0[i][j] = U1[i][j];
	  U1[i][j] = U2[i][j];
	}
      }

      // ピーク位置の追跡
      double max_temp = -INFINITY;
      int peak_x = -1;
      for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
	  if (U2[i][j] > max_temp) {
	    max_temp = U2[i][j];
	    peak_x = i;
	  }
        }
      }
      peak_positions[istep] = peak_x * dx;
      times[istep] = tnow;

      // 速度の計算
      if (istep > 0) {
        double displacement = peak_positions[istep] - peak_positions[istep - 1];
        double time_difference = times[istep] - times[istep - 1];
        velocities[istep - 1] = displacement / time_difference;
      }

      // tnowの更新
      tnow += dt;
      istep += 1;

      // ファイルへの速度の出力
      if (istep > 1) {
        FILE *velocity_file;
        velocity_file = fopen("velocities.txt", "w");
        for (int i = 0; i < istep - 1; ++i) {
	  fprintf(velocity_file, "%lf %lf\n", dx*(double)i, velocities[i]);
        }
        fclose(velocity_file);
      }


      if (istep%5==0){
	FILE* fp_phi;
	char name[120];
	sprintf(name,"hoge_phi-%03d.dat",istep);
	fp_phi = fopen(name,"w");
	for(i=0;i<N;i++){
	  for(j=0;j<N;j++) {
	    fprintf(fp_phi,"%12.4e %12.4e %12.4e\n",dx*(double)i,dx*(double)j,U2[i][j]);
	  }
	  fprintf(fp_phi,"\n");
	}
	fclose(fp_phi);
      }
    }

    for(i=0;i<N;i++){
      for(j=0;j<N;j++) {
	fprintf(output2_phi,"%12.4e %12.4e %12.4e\n",dx*(double)i,dx*(double)j,U2[i][j]);
      }
      fprintf(output2_phi,"\n");
    }
    fprintf(output2_phi,"\n");

    fclose(output2_phi);

    return 0;
  }
