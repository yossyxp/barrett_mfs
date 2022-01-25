
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//--------------------定数--------------------//

#define NUM ( 90 )
#define TMAX ( 600.0 )
//#define TMAX ( 1000 * dt )
//#define dt ( 0.1 / ( NUM * NUM ) ) // 0.000001
#define dt ( 1.0e-4 ) // 0.000001
#define RANGE_CHECK(x, xmin, xmax) ( x = ( x < xmin ? xmin : ( x < xmax ?  x : xmax)));

//----------定数(結晶)----------//

#define C ( 2.0e+7 )
#define rho ( 1.42e-3 )
#define alpha ( 1.0e-5 )
#define u_infty ( 0.004 )
#define eps ( 0.01 )
#define epsl ( 0.1 )
//#define theta0 ( M_PI / 12.0 )
#define theta0 ( 0.0 )

//--------------------関数--------------------//

double* make_vector( int N ); // 領域の確保(ベクトル)
double** make_matrix( int ROW, int COL ); // 領域の確保(行列)

void connect( double *x ); // 点のつなぎ(1個)
void connect_double( double* x, double *y ); // 点のつなぎ(2個)

double ip( double x1, double x2, double y1, double y2 ); // 内積
double DIST( double x1, double x2, double y1, double y2 ); // |y-x|
double DET( double x1, double x2, double y1, double y2 ); // 行列式

void gauss(int n, double **A, double *b, double *x); // Ax = b(0〜n-1までの行列)

void runge_kutta( double t, double *X1, double *X2 ); // ルンゲクッタ
void euler( double t, double *X1, double *X2 );

void F( double t, double *x1, double *x2, double *F1, double *F2 ); // 右辺
void pre( double t, double *x1, double *x2, double *T1, double *T2, double *N1, double *N2, double *V, double *W ); // 準備

void initial_condition( double *x1, double *x2 );

void quantities( double t, double *x1, double *x2, double *l, double *l_mid, double *t1, double *t2, double *n1, double *n2, double *T1, double *T2, double *N1, double *N2, double *nu, double *phi, double *kappa, double *kappa_mid, double *beta, double *kappa_gamma ); //x --> t,n,T,N,phi,kappa
void measure( double t, double *x1, double *x2, double *L, double *A ); // x,l --> L,A

double E( double x1, double x2, double y1, double y2 ); // 基本解
void grad_E( double x1, double x2, double y1, double y2, double *grad_E1, double *grad_E2 ); // 基本解の勾配

void PP( double t, double *beta, double *kappa_gamma, double *P, double A ); // P
void grad_P( double t, double *x_mid1, double *x_mid2, double *n1, double *n2, double *l, double A, double *P, double *beta, double *grad_u1, double *grad_u2 ); // x,n,r,P --> grad_aP

double omega( int n );

void velocity( double t, double *x1, double *x2, double *n1, double *n2, double *l, double *l_mid, double *phi, double *beta, double *kappa, double *kappa_mid, double *kappa_gamma, double L, double A, double *V, double *W ); // x,n,r,t,ohi,kappa,L --> V,W
void normal_speed( double t, double *x_mid1, double *x_mid2, double *n1, double *n2, double *phi, double *beta, double *kappa_gamma, double *grad_u1, double *grad_u2, double *v, double *V );
void tangent_speed( double t, double *l, double *phi, double *kappa, double *v, double *V, double L, double *W );
void tangent_speed_kappa( double t, double *l, double *l_mid, double *phi, double *kappa, double *kappa_mid, double *v, double *V, double L, double *W );

double phiphi( double x );
double phiphi_x( double x );
double Gamma( double l, int i, double *p1, double *p2 );
double Gamma_1( double l, int i, double *p1, double *p2 );
double Gamma_2( double l, int i, double *p1, double *p2 );
double Gamma_3( double l, int i );
double Gamma_4( double l, int i );
double Gammap1( double l, int i, double *p1, double *p2 );
double Gammap2( double l, int i, double *p1, double *p2 );
double Gamma_p1p1( double l, int i, double *p1, double *p2 );
double Gamma_p1p2( double l, int i, double *p1, double *p2 );
double Gamma_p2p2( double l, int i, double *p1, double *p2 );

//--------------------main-------------------//

int main(void){
  
  int i,z;
  double t;
  double *X1,*X2;
  X1 = make_vector(NUM + 2);
  X2 = make_vector(NUM + 2);
  
  double L,A;

  char file[5000];
  FILE *fp,*fp2;
  
  //fp = fopen("mfs.dat", "w");
  fp2 = fopen("L_A.dat", "w");
  
  t = 0.0;
  z = 0;

  initial_condition(X1,X2);
  
  sprintf(file, "./data/barrett_mfs%06d.dat", z);
  fp = fopen(file, "w");
  
  for( i = 0; i <= NUM; i++ ){
    
    fprintf(fp, "%f %f %f\n", X1[i], X2[i], t);
    
  }
  fprintf(fp,"\n");
  
  measure(t,X1,X2,&L,&A);
  
  fprintf(fp2, "t = %f, L = %f, A = %f\n", t, L, A );

  fclose(fp2);

  fclose(fp);
  
  while( t < TMAX ){
    
    runge_kutta(t,X1,X2);
    //euler(t,X1,X2);
    
    z++;
    t = z * dt;

    measure(t,X1,X2,&L,&A);
    printf("t = %f A = %.15f\n", t, A);
    
    if( z % 100 == 0 ){
      
      sprintf(file, "./data/barrett_mfs%06d.dat", z / 100 );
      fp = fopen(file, "w");

      measure(t,X1,X2,&L,&A);
      
      for( i = 0; i <= NUM; i++ ){
	
	fprintf(fp, "%f %f %f\n", X1[i], X2[i], t);
	
      }
      fprintf(fp,"\n");
      
      fclose(fp);
      fp2 = fopen("L_A.dat", "a");
      
      fprintf(fp2, "t = %f, L = %f, A = %f\n", t, L, A );
      fclose(fp2);
      
    }

  }
  
  free(X1);
  free(X2);
  
  return 0;
  
}


//--------------------関数--------------------//

// 領域の確保(ベクトル)
double* make_vector( int N ){
  
  double *a;
  
  // メモリの確保
  if( ( a = malloc( sizeof( double ) * N ) ) == NULL ){
    
    printf( "LACK of AVAILABLE MEMORY!" );
    exit(1);

  }
  
  return a;

}

// 領域の確保(行列)
double** make_matrix(int ROW, int COL){
  
  int i;
  double **b;
  
  // メモリの確保
  if( ( b = malloc( sizeof( double* ) * ROW ) ) == NULL ){
    
    printf( "LACK of AVAILABLE MEMORY!" );
    exit(1);

  }
  
  for( i = 0; i < ROW; i++){
    
    if( ( b[i] = malloc( sizeof( double ) * COL ) ) == NULL ){
      
      printf( "LACK of AVAILABLE MEMORY!" );
      free(b);
      exit(1);

    }

  }
  
  return b;

}

// 点のつなぎ(1個)
void connect( double *x ){
  
  x[0] = x[NUM];
  x[NUM + 1] = x[1];

}

// 点のつなぎ(2個)
void connect_double( double *x, double *y ){
  
  connect(x);
  connect(y);

}

// 内積
double ip( double x1, double x2, double y1, double y2 ){
  
  return ( x1 * y1 + x2 * y2 );

}

// |y-x|
double DIST(double x1, double x2, double y1, double y2){
  
  return ( sqrt( (y1 - x1) * (y1 - x1) + (y2 - x2) * (y2 - x2) ) );

}

// 行列式
double DET( double x1, double x2, double y1, double y2 ){
  
  return ( x1 * y2 - y1 * x2 );
  
}

// Ax = b(0〜n-1までの行列)
void gauss( int n, double **A, double *b, double *x ){
  
  int i,j,k,l,row_max;
  double max,temp,c;
  double temp_vec;
  double *temp_mat;
  
  for( i = 0; i < n - 1; i++ ){
    
    row_max = i;
    max = A[i][i];
    
    for( l = i + 1; l < n; l++ ){
      
      if( max < A[l][i] ){
	
	row_max = l;
	max = A[l][i];

      }		    

    }
    
    if( row_max > i ){
      
      temp_mat = A[i];
      A[i] = A[row_max];
      A[row_max] = temp_mat;
      
      temp_vec = b[i];
      b[i] = b[row_max];
      b[row_max] = temp_vec;

    }
      
    for( k = i + 1; k < n; k++ ){
      
      c = A[k][i] / A[i][i];

      for(j = i; j < n; j++){
	
	A[k][j] -= A[i][j] * c;

      }

      b[k] -= b[i]*c;

    }

  }
  
  for( i = n - 1; i >= 0; i--){
    
    temp = b[i];
    
    for( j = n - 1; j > i; j-- ){
      
      temp -= A[i][j] * x[j];

    }

    x[i] = temp / A[i][i];
    
  }

}

// ルンゲクッタ
void runge_kutta( double t, double *X1, double *X2 ){
  
  int i;
  double t_temp;
  double *x_temp1,*x_temp2,*F1,*F2;
  double *k11,*k12,*k21,*k22,*k31,*k32,*k41,*k42;
  
  k11 = make_vector(NUM + 2);
  k12 = make_vector(NUM + 2);
  k21 = make_vector(NUM + 2);
  k22 = make_vector(NUM + 2);
  k31 = make_vector(NUM + 2);
  k32 = make_vector(NUM + 2);
  k41 = make_vector(NUM + 2);
  k42 = make_vector(NUM + 2);

  x_temp1 = make_vector(NUM + 2);
  x_temp2 = make_vector(NUM + 2);
  F1 = make_vector(NUM + 2);
  F2 = make_vector(NUM + 2);

  F(t,X1,X2,F1,F2);
  
  for( i = 1; i <= NUM; i++ ){
    
      k11[i] = F1[i];
      k12[i] = F2[i];
      x_temp1[i] = X1[i] + dt * k11[i] / 2.0;
      x_temp2[i] = X2[i] + dt * k12[i] / 2.0;

  }
  connect_double(x_temp1,x_temp2);

  t_temp = t + dt/2.0;

  F(t_temp,x_temp1,x_temp2,F1,F2);      

  for( i = 1; i <= NUM; i++ ){
    
      k21[i] = F1[i];
      k22[i] = F2[i];
      x_temp1[i] = X1[i] + dt * k21[i] / 2.0;
      x_temp2[i] = X2[i] + dt * k22[i] / 2.0;

  }
  connect_double(x_temp1,x_temp2);

  F(t_temp,x_temp1,x_temp2,F1,F2);
  
  for( i = 1; i <= NUM; i++ ){
    
    k31[i] = F1[i];
    k32[i] = F2[i];
    x_temp1[i] = X1[i] + k31[i] * dt;
    x_temp2[i] = X2[i] + k32[i] * dt;

  }
  connect_double(x_temp1,x_temp2);
  
  t_temp = t + dt;
      
  F(t_temp,x_temp1,x_temp2,F1,F2);
  
  for( i = 1; i <= NUM; i++ ){
    
    k41[i] = F1[i];
    k42[i] = F2[i];
    
    X1[i] = X1[i] + dt * ( k11[i] + 2.0 * k21[i] + 2.0 * k31[i] + k41[i] ) / 6.0;
    X2[i] = X2[i] + dt * ( k12[i] + 2.0 * k22[i] + 2.0 * k32[i] + k42[i] ) / 6.0;

  }
  connect_double(X1,X2);

  free(k11);
  free(k12);
  free(k21);
  free(k22);
  free(k31);
  free(k32);
  free(k41);
  free(k42);
  free(x_temp1);
  free(x_temp2);
  free(F1);
  free(F2);
  
}

void euler( double t, double *X1, double *X2 ){

  int i;
  double *F1,*F2;

  F1 = make_vector(NUM + 2);
  F2 = make_vector(NUM + 2);

  F(t,X1,X2,F1,F2);

  for( i = 1; i <= NUM; i++ ){
    
    X1[i] = X1[i] + dt * F1[i];
    X2[i] = X2[i] + dt * F2[i];

  }
  connect_double(X1,X2);
  
  free(F1);
  free(F2);
  
}

// 右辺
void F( double t, double *x1, double *x2, double *F1, double *F2 ){
  
  int i;
  
  double *T1,*T2,*N1,*N2;
  T1 = make_vector(NUM + 2);
  T2 = make_vector(NUM + 2);
  N1 = make_vector(NUM + 2);
  N2 = make_vector(NUM + 2);

  double *V,*W;
  V = make_vector(NUM + 2);
  W = make_vector(NUM + 2);
  
  pre(t,x1,x2,T1,T2,N1,N2,V,W);
  
  for( i = 1; i <= NUM; i++ ){
    
    F1[i] = V[i] * N1[i] + W[i] * T1[i];
    F2[i] = V[i] * N2[i] + W[i] * T2[i];
    //printf("V = %f\n", V[i]);

  }
  connect_double(F1,F2);

  free(W);
  free(V);
  free(T1);
  free(T2);
  free(N1);
  free(N2);
  
}

//x --> T,N,V,W
void pre( double t, double *x1, double *x2, double *T1, double *T2, double *N1, double *N2, double *V, double *W ){
  
  double *l;
  double *l_mid;
  double *t1,*t2,*n1,*n2;
  double *nu;
  double *phi;
  double *beta;
  double *kappa;
  double *kappa_mid;
  double *kappa_gamma;
  double L,A;

  l = make_vector(NUM + 2);
  l_mid = make_vector(NUM + 2);
  nu = make_vector(NUM + 2);
  phi = make_vector(NUM + 2);
  kappa = make_vector(NUM + 2);
  kappa_mid = make_vector(NUM + 2);
  beta = make_vector(NUM + 2);
  kappa_gamma = make_vector(NUM + 2);

  t1 = make_vector(NUM + 2);
  t2 = make_vector(NUM + 2);
  n1 = make_vector(NUM + 2);
  n2 = make_vector(NUM + 2);
  
  // T,N
  quantities(t,x1,x2,l,l_mid,t1,t2,n1,n2,T1,T2,N1,N2,nu,phi,kappa,kappa_mid,beta,kappa_gamma);
  
  // L,A
  measure(t,x1,x2,&L,&A);
  
  // V,W
  velocity(t,x1,x2,n1,n2,l,l_mid,phi,beta,kappa,kappa_mid,kappa_gamma,L,A,V,W);
  
  free(t1);
  free(t2);
  free(n1);
  free(n2);

  free(kappa_gamma);
  free(beta);
  free(kappa_mid);
  free(kappa);
  free(phi);
  free(nu);
  free(l_mid);
  free(l);
  
}

// 初期曲線
void initial_condition( double *x1, double *x2 ){

  int i,k;
  double u;
  double lambda;

  
  for( i = 1; i <= NUM; i++ ){
    
    u = i * 2 * M_PI / NUM;
    
    x1[i] = 5.0e-2 * cos(u);
    x2[i] = 5.0e-2 * sin(u);
    
  }
  connect_double(x1,x2);
  
  /*
  for (i = 1; i <= NUM; i++)
  {

    u = i * 2 * M_PI / NUM;
  }

  for (k = 0; k < 6; k++)
  {
    for (i = k * (NUM / 6) + 1; i <= (k + 1) * (NUM / 6); i++)
    {
      lambda = (i - k * (NUM / 6)) * 1.0 / (NUM / 6);
      if (k == 0)
      {
        x1[i] = (1.0 - lambda) * 1.0 + lambda * 1.0 / 2.0;
        x2[i] = (1.0 - lambda) * 0.0 + lambda * sqrt(3.0) / 2.0;
      }
      else if (k == 1)
      {
        x1[i] = (1.0 - lambda) * 1.0 / 2.0 + lambda * (-1.0 / 2.0);
        x2[i] = sqrt(3.0) / 2.0;
      }
      else if (k == 2)
      {
        x1[i] = (1.0 - lambda) * (-1.0 / 2.0) + lambda * (-1.0);
        x2[i] = (1.0 - lambda) * (sqrt(3.0) / 2.0) + lambda * 0.0;
      }
      else if (k == 3)
      {
        x1[i] = (1.0 - lambda) * (-1.0) + lambda * (-1.0 / 2.0);
        x2[i] = (1.0 - lambda) * 0.0 + lambda * (-sqrt(3.0) / 2.0);
      }
      else if (k == 4)
      {
        x1[i] = (1.0 - lambda) * (-1.0 / 2.0) + lambda * (1.0 / 2.0);
        x2[i] = -sqrt(3.0) / 2.0;
      }
      else
      {
        x1[i] = (1.0 - lambda) * (1.0 / 2.0) + lambda * (1.0);
        x2[i] = (1.0 - lambda) * (-sqrt(3.0) / 2.0) + lambda * 0.0;
      }

      x1[i] = 5.0e-2 * x1[i];
      x2[i] = 5.0e-2 * x2[i];
    }
  }
  connect_double(x1, x2);
  */
  
}

//x --> t,n,T,N,phi,kappa
void quantities( double t, double *x1, double *x2, double *l, double *l_mid, double *t1, double *t2, double *n1, double *n2, double *T1, double *T2, double *N1, double *N2, double *nu, double *phi, double *kappa, double *kappa_mid, double *beta, double *kappa_gamma ){
  
  int i;
  double D,I,D_sgn;
  
  for( i = 1; i <= NUM; i++ ){
    
    l[i] = DIST(x1[i - 1],x2[i - 1],x1[i],x2[i]);
    
    t1[i] = ( x1[i] - x1[i - 1] ) / l[i];
    t2[i] = ( x2[i] - x2[i - 1] ) / l[i];
    
    n1[i] = t2[i];
    n2[i] = -t1[i];

  }
  connect(l);
  connect_double(t1,t2);
  connect_double(n1,n2);
  

  RANGE_CHECK(t1[1],-1.0,1.0);
  
  if( t2[1] >= 0 ){
    
    nu[1] = acos(t1[1]);
    
  }
  
  else{
    
    nu[1] = -acos(t1[1]);
    
  }
  
  for( i = 1; i <= NUM; i++){
    
    D = DET(t1[i],t2[i],t1[i + 1],t2[i + 1]);
    I = ip(t1[i],t2[i],t1[i + 1],t2[i + 1]);
    
    RANGE_CHECK(I,-1.0,1.0);
    
    if( D < 0 ){
      
      D_sgn = -1;
      
    }
    
    else if( D > 0 ){
      
      D_sgn = 1;
      
    }
    
    else{
      
      D_sgn = 0;
      
    }
    
    nu[i + 1] = nu[i] + D_sgn * acos(I);
    
  }
  nu[0] = nu[1] - ( nu[NUM + 1] - nu[NUM] );
  
  
  for( i = 1; i <= NUM; i++ ){
    
    phi[i] = nu[i + 1] - nu[i];
    l_mid[i] = ( l[i] + l[i + 1] ) / 2.0;
    
  }
  connect(phi);
  connect(l_mid);
  

  for( i = 1; i <= NUM; i++){
    
    T1[i] = ( t1[i] + t1[i + 1] ) / ( 2.0 * cos(phi[i] / 2.0) );
    T2[i] = ( t2[i] + t2[i + 1] ) / ( 2.0 * cos(phi[i] / 2.0) );

    N1[i] = T2[i];
    N2[i] = -T1[i];
    
    kappa[i] = ( tan(phi[i] / 2.0) + tan(phi[i - 1] / 2.0) ) / l[i];

    kappa_mid[i] = 2 * sin(phi[i] / 2.0) / l_mid[i];
    
  }
  connect_double(T1,T2);
  connect_double(N1,N2);
  connect(kappa);
  connect(kappa_mid);
  
  for( i = 1; i <= NUM; i++ ){
    
    //beta[i] = 1.0;
    //beta[i] = ( 1 + eps * cos(6 * nu[i]) );

    /*
    beta[i] = sqrt(( cos(M_PI/3.0)*n1[i] + sin(M_PI/3.0)*n2[i] ) * ( cos(M_PI/3.0)*n1[i] + sin(M_PI/3.0)*n2[i] ) + eps * eps * ( -sin(M_PI/3.0)*n1[i] + cos(M_PI/3.0)*n2[i] ) * ( -sin(M_PI/3.0)*n1[i] + cos(M_PI/3.0)*n2[i] ))
            + sqrt(( cos(2*M_PI/3.0)*n1[i] + sin(2*M_PI/3.0)*n2[i] ) * ( cos(2*M_PI/3.0)*n1[i] + sin(2*M_PI/3.0)*n2[i] ) + eps * eps * ( -sin(2*M_PI/3.0)*n1[i] + cos(2*M_PI/3.0)*n2[i] ) * ( -sin(2*M_PI/3.0)*n1[i] + cos(2*M_PI/3.0)*n2[i] ))
            + sqrt(( cos(M_PI)*n1[i] + sin(M_PI)*n2[i] ) * ( cos(M_PI)*n1[i] + sin(M_PI)*n2[i] ) + eps * eps * ( -sin(M_PI)*n1[i] + cos(M_PI)*n2[i] ) * ( -sin(M_PI)*n1[i] + cos(M_PI)*n2[i] ));
    */
    
    beta[i] = Gamma(1.0,i,n1,n2) + Gamma(2.0,i,n1,n2) + Gamma(3.0,i,n1,n2);
    
    //printf("%.15f\n", beta[i]);
  }
  connect(beta);

  for( i = 1; i <= NUM; i++ ){

    //kappa_gamma[i] = kappa[i];
    //kappa_gamma[i] = -( beta[i] - ( 36 * eps * cos(6 * nu[i]) ) ) * kappa[i];
    //kappa_gamma[i] = ( beta[i] + ( 0.99 / 35.0 ) * 36 * cos(6 * nu[i]) );
    
    
    //kappa_gamma[i] = -( Gamma_p1p1(1.0,i,n1,n2) + Gamma_p1p1(2.0,i,n1,n2) + Gamma_p1p1(3.0,i,n1,n2) + Gamma_p2p2(1.0,i,n1,n2) + Gamma_p2p2(2.0,i,n1,n2) + Gamma_p2p2(3.0,i,n1,n2) );

    //kappa_gamma[i] = -( ( Gamma_p1p1(1.0,i,n1,n2) + Gamma_p1p1(2.0,i,n1,n2) + Gamma_p1p1(3.0,i,n1,n2) ) * n2[i] * n2[i] + 2 * ( Gamma_p1p2(1.0,i,n1,n2) + Gamma_p1p2(2.0,i,n1,n2) + Gamma_p1p2(3.0,i,n1,n2) ) * n1[i] * n2[i] + ( Gamma_p2p2(1.0,i,n1,n2) + Gamma_p2p2(2.0,i,n1,n2) + Gamma_p2p2(3.0,i,n1,n2) ) * n1[i] * n1[i] ) * kappa[i];
    
    kappa_gamma[i] = -( beta[i] + ( Gamma_p1p1(1.0,i,n1,n2) + Gamma_p1p1(2.0,i,n1,n2) + Gamma_p1p1(3.0,i,n1,n2) ) * cos(nu[i]) * cos(nu[i]) + 2 * ( Gamma_p1p2(1.0,i,n1,n2) + Gamma_p1p2(2.0,i,n1,n2) + Gamma_p1p2(3.0,i,n1,n2) ) * cos(nu[i]) * sin(nu[i]) + ( Gamma_p2p2(1.0,i,n1,n2) + Gamma_p2p2(2.0,i,n1,n2) + Gamma_p2p2(3.0,i,n1,n2) ) * sin(nu[i]) * sin(nu[i]) + ( Gammap1(1.0,i,n1,n2) + Gammap1(2.0,i,n1,n2) + Gammap1(3.0,i,n1,n2) ) * ( -sin(nu[i]) ) + ( Gammap2(1.0,i,n1,n2) + Gammap2(2.0,i,n1,n2) + Gammap2(3.0,i,n1,n2) ) * cos(nu[i]) ) * kappa[i];
    
    /*
    kappa_gamma[i] = -( beta[i] + ( ( Gamma_p1p1(1.0,i,n1,n2) + Gamma_p1p1(2.0,i,n1,n2) + Gamma_p1p1(3.0,i,n1,n2) ) * cos(nu[i]) + ( Gamma_p1p2(1.0,i,n1,n2) + Gamma_p1p2(2.0,i,n1,n2) + Gamma_p1p2(3.0,i,n1,n2) ) * sin(nu[i]) ) * cos(nu[i]) + ( Gammap1(1.0,i,n1,n2) + Gammap1(2.0,i,n1,n2) + Gammap1(3.0,i,n1,n2) ) * ( -sin(nu[i]) ) + ( ( Gamma_p1p2(1.0,i,n1,n2) + Gamma_p1p2(2.0,i,n1,n2) + Gamma_p1p2(3.0,i,n1,n2) ) * cos(nu[i]) + ( Gamma_p2p2(1.0,i,n1,n2) + Gamma_p2p2(2.0,i,n1,n2) + Gamma_p2p2(3.0,i,n1,n2) ) * sin(nu[i]) ) * sin(nu[i]) + ( Gammap2(1.0,i,n1,n2) + Gammap2(2.0,i,n1,n2) + Gammap2(3.0,i,n1,n2) ) * cos(nu[i]) ) * kappa[i];
    */
    
    //printf("%f %f %.15f\n", n1[i], n2[i], kappa_gamma[i]);
    
  }
  connect(kappa_gamma);
  
}

// x,l --> L,A
void measure( double t, double *x1, double *x2, double *L, double *A ){
  
  int i;
  
  *L = 0.0;
  *A = 0.0;
  
  for( i = 1; i <= NUM; i++){
    
    *L += DIST(x1[i],x2[i],x1[i - 1],x2[i - 1]);
    *A += DET(x1[i - 1],x2[i - 1],x1[i],x2[i]);

  }
  *A = *A / 2.0;
  
}

// 基本解
double E( double x1, double x2, double y1, double y2 ){
  
  return ( log(DIST(y1,y2,x1,x2)) / ( 2.0 * M_PI ) );

}

// 基本解の勾配
void grad_E( double x1, double x2, double y1, double y2, double *grad_E1, double *grad_E2 ){
  
  double ry = DIST(y1,y2,x1,x2);
  
  *grad_E1 = ( x1 - y1 ) / ( 2.0 * M_PI * ry * ry );
  *grad_E2 = ( x2 - y2 ) / ( 2.0 * M_PI * ry * ry );
  
}

// P
void PP( double t, double *beta, double *kappa_gamma, double *P, double A ){
  
  int i;
  double r_c,R;
  
  r_c = sqrt(A / M_PI);
  R = 6.5 * r_c;

  P[0] = 0.0;
  //P[0] = 2 * M_PI * sigma_infty / log(R / r_c);
  
  for( i = 1; i <= NUM; i++ ){
    
    P[i] = -alpha * kappa_gamma[i];
    
  }

  for( i = NUM + 1; i <= 2 * NUM; i++ ){
    
    P[i] = u_infty;
    
  }
  
}

// x,n,r,P --> grad_P
void grad_P( double t, double *x_mid1, double *x_mid2, double *n1, double *n2, double *l, double A, double *P, double *beta, double *grad_u1, double *grad_u2 ){
  
  int i,j,k;
  double r_c,R;
  double lambda;
  double theta;
  double amano;
  double *x3,*x4;
  double *ll;
  double *tt1,*tt2;
  double *nn1,*nn2;
  double *x_mid3,*x_mid4;
  double *y1,*y2,*z1,*z2;
  double **G;
  double **H1,**H2;
  double *H;
  double *Q;
  double grad_E1,grad_E2;
  double *amano1,*amano2,*amano3,*amano4;
  
  //double d = 1.0 / sqrt(1.0 * NUM);
  
  double d;
  
  x3 = make_vector(NUM + 2);
  x4 = make_vector(NUM + 2);
  x_mid3 = make_vector(NUM + 2);
  x_mid4 = make_vector(NUM + 2);
  y1 = make_vector(2 * NUM + 2);
  y2 = make_vector(2 * NUM + 2);
  z1 = make_vector(2 * NUM + 2);
  z2 = make_vector(2 * NUM + 2);
  
  amano1 = make_vector(NUM + 2);
  amano2 = make_vector(NUM + 2);
  amano3 = make_vector(NUM + 2);
  amano4 = make_vector(NUM + 2);
  
  ll = make_vector(NUM + 2);
  tt1 = make_vector(NUM + 2);
  tt2 = make_vector(NUM + 2);
  nn1 = make_vector(NUM + 2);
  nn2 = make_vector(NUM + 2);
  
  H1 = make_matrix(2 * NUM + 1,2 * NUM + 1);
  H2 = make_matrix(2 * NUM + 1,2 * NUM + 1);

  H = make_vector(2 * NUM + 1);
  Q = make_vector(2 * NUM + 1);

  G = make_matrix(2 * NUM + 1,2 * NUM + 1);

  
  r_c = sqrt(A / M_PI);
  R = 6.5 * r_c;

  
  for( i = 1; i <= NUM; i++ ){
    
    theta = i * 2 * M_PI / NUM;
    
    x3[i] = 4.0 * cos(theta);
    x4[i] = 4.0 * sin(theta);

    //printf("%f %f\n", x3[i], x4[i]);
    
  }
  connect_double(x3,x4);
  
  
  /*
  for( k = 0; k < 8; k++ ){
    
    for( i = k * ( NUM / 8 ) + 1; i <= ( k + 1 ) * ( NUM / 8 ); i++ ){
      
      lambda = ( i - k * ( NUM / 8 ) ) * 1.0 / ( NUM / 8 );

      if( k == 0 ){
	
        x3[i] = 1.0;
        x4[i] = lambda;
	
      }
      
      else if( k == 1 ){
	
        x3[i] = 1.0 - lambda;
        x4[i] = 1.0;
	
      }
      
      else if( k == 2 ){
	
	x3[i] = 0.0 - lambda;
        x4[i] = 1.0;

      }
      
      else if( k == 3 ){

	x3[i] = -1.0;
        x4[i] = 1.0 - lambda;

      }

      else if( k == 4 ){
	
	x3[i] = -1.0;
        x4[i] = 0.0 - lambda;

      }

      else if( k == 5 ){
	
	x3[i] = -1.0 + lambda;
        x4[i] = -1.0;

      }

      else if( k == 6 ){
	
	x3[i] = 0.0 + lambda;
        x4[i] = -1.0;

      }

      else if( k == 7 ){
	
	x3[i] = 1.0;
        x4[i] = -1.0 + lambda;

      }

      x3[i] = 4.0 * x3[i];
      x4[i] = 4.0 * x4[i];
      
      //printf("%f %f\n", x3[i], x4[i]);
      
    }
    
  }
  connect_double(x3, x4);
  */
  
  
  for( i = 1; i <= NUM; i++ ){
    
    ll[i] = DIST(x3[i - 1],x4[i - 1],x3[i],x4[i]);
    
    tt1[i] = ( x3[i] - x3[i - 1] ) / ll[i];
    tt2[i] = ( x4[i] - x4[i - 1] ) / ll[i];
    
    nn1[i] = tt2[i];
    nn2[i] = -tt1[i];

  }
  connect(ll);
  connect_double(tt1,tt2);
  connect_double(nn1,nn2);  
  

  for( i = 1; i <= NUM; i++ ){
    
    x_mid3[i] = ( x3[i] + x3[i - 1] ) / 2.0;
    x_mid4[i] = ( x4[i] + x4[i - 1] ) / 2.0;
    
  }
  connect_double(x_mid3,x_mid4);

  for( i = 1; i <= NUM; i++ ){
    
    amano1[i] = -( -( x_mid2[i + 1] - x_mid2[i - 1]) / DIST(x_mid1[i - 1],x_mid2[i - 1],x_mid1[i + 1],x_mid2[i + 1]) );
    
    amano2[i] = - ( ( x_mid1[i + 1] - x_mid1[i - 1]) / DIST(x_mid1[i - 1],x_mid2[i - 1],x_mid1[i + 1],x_mid2[i + 1]) );
    
    amano3[i] = -( -( x_mid4[i + 1] - x_mid4[i - 1]) / DIST(x_mid3[i - 1],x_mid4[i - 1],x_mid3[i + 1],x_mid4[i + 1]) );
    
    amano4[i] = -( ( x_mid3[i + 1] - x_mid3[i - 1]) / DIST(x_mid3[i - 1],x_mid4[i - 1],x_mid3[i + 1],x_mid4[i + 1]) );

    //printf("%f %f\n", amano1[i], amano2[i]);
    
  }
  connect_double(amano1,amano2);
  connect_double(amano3,amano4);
  
  for( i = 1; i <= NUM; i++ ){

    d = DIST(x_mid1[i - 1],x_mid2[i - 1],x_mid1[i + 1],x_mid2[i + 1]) / 2.0;
    
    //y1[i] = x_mid1[i] - d * n1[i];
    //y2[i] = x_mid2[i] - d * n2[i];

    y1[i] = x_mid1[i] - d * amano1[i]; 
    y2[i] = x_mid2[i] - d * amano2[i];

    //z1[i] = 0.0;
    //z2[i] = 0.0;

    //printf("%f %f %d\n", y1[i], y2[i], i);
    
  }
  
  for( i = NUM + 1; i <= 2 * NUM; i++ ){

    d = DIST(x_mid3[i - 1 - NUM],x_mid4[i - 1 - NUM],x_mid3[i + 1 - NUM],x_mid4[i + 1 - NUM]) / 2.0;

    //y1[i] = x_mid3[i - NUM] + d * nn1[i - NUM];
    // y2[i] = x_mid4[i - NUM] + d * nn2[i - NUM];
    
    y1[i] = x_mid3[i - NUM] + d * amano3[i - NUM]; 
    y2[i] = x_mid4[i - NUM] + d * amano4[i - NUM];

    //z1[i] = 1000.0 * y1[i];
    //z2[i] = 1000.0 * y2[i];

    //printf("%f %f\n", x_mid3[i - NUM], x_mid4[i - NUM]);
    //printf("%f %f\n", y1[i], y2[i]);    
    
  }

  for( i = 1; i <= NUM; i++ ){
    
    for( j = 1; j <= 2 * NUM; j++ ){
      
      grad_E(x_mid1[i],x_mid2[i],y1[j],y2[j],&grad_E1,&grad_E2);
      
      H1[i][j] = grad_E1;
      H2[i][j] = grad_E2;
      
    }
    
  }

  for( i = NUM + 1; i <= 2 * NUM; i++ ){
    
    for( j = 1; j <= 2 * NUM; j++ ){
      
      grad_E(x_mid3[i - NUM],x_mid4[i - NUM],y1[j],y2[j],&grad_E1,&grad_E2);
      
      H1[i][j] = grad_E1;
      H2[i][j] = grad_E2;
      
    }

  }

  /*
  for( j = 1; j <= 2 * NUM; j++){
    
    H[j] = 0.0;
    
    for( i = 1; i <= NUM; i++){
      
      H[j] += ip(H1[i][j],H2[i][j],n1[i],n2[i]) * l[i];

    }
    
  }
  */
  
  G[0][0] = 0.0;

  for( j = 1; j <= 2 * NUM; j++ ){
    
    G[0][j] = 1.0;
    
  }

  for( i = 1; i <= 2 * NUM; i++ ){
    
    G[i][0] = 1.0;
    
  }
  
  for( i = 1; i <= NUM; i++ ){
    
    for( j = 1; j <= 2 * NUM; j++ ){
      
      G[i][j] = E(x_mid1[i],x_mid2[i],y1[j],y2[j]) - rho * ip(H1[i][j],H2[i][j],n1[i],n2[i]) / beta[i];
      
    }
    
  }
  
  for( i = NUM + 1; i <= 2 * NUM; i++ ){
    
    for( j = 1; j <= 2 * NUM; j++ ){
      
      G[i][j] = E(x_mid3[i - NUM],x_mid4[i - NUM],y1[j],y2[j]);
      
    }
    
  }

  /*
  for( i = 0; i <= 2 * NUM; i++ ){
    
    for( j = 0; j <= 2 * NUM; j++ ){
    
    printf("G[%d][%d] = %.30f\n", i, j, G[i][j]);

    }
    
  }
  */
  
  gauss(2 * NUM + 1,G,P,Q);

  /*
  for( i = 1; i <= NUM; i++ ){

    printf("u[%d] = %.15f\n", i, Q[i]);

  }
  */
  
  for( i = 1; i <= NUM; i++ ){

    grad_u1[i] = 0.0;
    grad_u2[i] = 0.0;

    for( j = 1; j <= 2 * NUM; j++ ){
      
      grad_u1[i] += Q[j] * H1[i][j];
      grad_u2[i] += Q[j] * H2[i][j];
      
    }
    
    //printf("grad_u1[%d] = %.15f grad_u2[%d] = %.15f\n", i, grad_u1[i], i, grad_u2[i]);
    
  }
  connect_double(grad_u1,grad_u2);

  /*
  for( i = 1; i <= NUM; i++ ){
    
    printf("u[%d] = %.15f\n", i, u[i]);
    
  }
  */
  
  for( i = 0; i <= 2 * NUM; i++ ){
    
    free(G[i]);
    
  }
  free(G);

  free(x3);
  free(x4);
  free(x_mid3);
  free(x_mid4);
  free(y1);
  free(y2);
  free(z1);
  free(z2);
  
  free(amano1);
  free(amano2);
  free(amano3);
  free(amano4);

  free(ll);
  free(tt1);
  free(tt2);
  free(nn1);
  free(nn2);
  
  for( i = 0; i <= 2 * NUM; i++ ){
    
    free(H1[i]);
    
  }
  free(H1);
  
  for( i = 0; i <= 2 * NUM; i++ ){
    
    free(H2[i]);
    
  }
  free(H2);

  free(Q);
  free(H);
  
}

double omega( int n ){
  
  return ( 10.0 * n );
  //return ( 0.1 / dt );

}

void velocity( double t, double *x1, double *x2, double *n1, double *n2, double *l, double *l_mid, double *phi, double *beta, double *kappa, double *kappa_mid, double *kappa_gamma, double L, double A, double *V, double *W  ){
  
  int i;
  double *P;
  double *x_mid1,*x_mid2;
  double *v;
  double grad_E1, grad_E2;
  double *grad_u1, *grad_u2;
  
  x_mid1 = make_vector(NUM + 2);
  x_mid2 = make_vector(NUM + 2);
  grad_u1 = make_vector(NUM + 2);
  grad_u2 = make_vector(NUM + 2);
  P = make_vector(2 * NUM + 1);
  v = make_vector(NUM + 2);
  
  for( i = 1; i <= NUM; i++ ){
    
    x_mid1[i] = ( x1[i] + x1[i - 1] ) / 2.0;
    x_mid2[i] = ( x2[i] + x2[i - 1] ) / 2.0;

  }
  connect_double(x_mid1,x_mid2);
  
  PP(t,beta,kappa_gamma,P,A);
  grad_P(t,x_mid1,x_mid2,n1,n2,l,A,P,beta,grad_u1,grad_u2);
  normal_speed(t,x_mid1,x_mid2,n1,n2,phi,beta,kappa_gamma,grad_u1,grad_u2,v,V);
  tangent_speed(t,l,phi,kappa,v,V,L,W);
  //tangent_speed_kappa(t,l,l_mid,phi,kappa,kappa_mid,v,V,L,W);

  free(x_mid1);
  free(x_mid2);
  free(grad_u1);
  free(grad_u2);
  free(v);
  free(P);
  
}

void normal_speed( double t, double *x_mid1, double *x_mid2, double *n1, double *n2, double *phi, double *beta, double *kappa_gamma, double *grad_u1, double *grad_u2, double *v, double *V ){
  
  int i;
  
  for( i = 1; i <= NUM; i++ ){

    //v[i] = beta[i];
    v[i] = ip(grad_u1[i],grad_u2[i],n1[i],n2[i]);

    //printf("v[%d] = %.15f\n", i, v[i]);

  }
  connect(v);

  for( i = 1; i <= NUM; i++ ){
    
    V[i] = ( v[i] + v[i + 1] ) / ( 2.0 * cos(phi[i] / 2.0) );
  }
  connect(V);
  
}

void tangent_speed( double t, double *l, double *phi, double *kappa, double *v, double *V, double L, double *W ){
  
  int i;
  double *psi,*PSI;
  double L_dot;
  double a,b,c;
  
  psi = make_vector(NUM + 1);
  PSI = make_vector(NUM + 1);
  
  psi[1] = 0.0;
  L_dot = 0.0;
  
  for(i = 1; i <= NUM; i++ ){
    
    L_dot += kappa[i] * v[i] * l[i];

  }
  
  for( i = 2; i <= NUM; i++ ){
    
    psi[i] = ( L_dot / NUM ) - V[i] * sin(phi[i] / 2.0) - V[i-1] * sin(phi[i - 1] / 2.0) + ( ( L / NUM ) - l[i] ) * omega(NUM);

  }
  
  PSI[1] = psi[1];
  
  for( i = 2; i <= NUM; i++ ){
    
    PSI[i] = PSI[i-1] + psi[i];

  }

  a = 0.0;
  b = 0.0;

  for( i = 1; i <= NUM; i++ ){
    
    a += PSI[i] / cos(phi[i] / 2.0);
    b += 1.0 / cos(phi[i] / 2.0);

  }
  c = -a / b;

  for( i = 1; i <= NUM; i++ ){
    
    W[i] = ( PSI[i] + c ) / cos(phi[i] / 2.0);

  }
  connect(W);

  free(PSI);
  free(psi);

}

void tangent_speed_kappa( double t, double *l, double *l_mid, double *phi, double *kappa, double *kappa_mid, double *v, double *V, double L, double *W ){
  
  int i,j;
  double *psi,*PSI;
  double *v_x;
  double *f;
  double f_kakko,phi_kakko,V_kakko;
  double L_dot;
  double a,b,c;
  
  psi = make_vector(NUM + 1);
  PSI = make_vector(NUM + 1);
  v_x = make_vector(NUM + 2);
  f = make_vector(NUM + 2);

  for(i = 1; i <= NUM; i++ ){

    v_x[i] = ( v[i + 1] - v[i] ) / l_mid[i];
    
  }
  connect(v_x);
  
  for(i = 1; i <= NUM; i++ ){

    f[i] = phiphi(kappa[i]) * kappa[i] * v[i] - phiphi_x(kappa[i]) * ( ( ( v_x[i] - v_x[i - 1] ) / l[i] ) + kappa[i] * kappa[i] * v[i] );

  }
  connect(f);

  f_kakko = 0.0;
  for( i = 1; i <= NUM; i++ ){

    f_kakko += f[i] * l[i];
    
  }
  f_kakko = f_kakko / L;

  phi_kakko = 0.0;
  for( i = 1; i <= NUM; i++ ){

    phi_kakko += phiphi(kappa[i]) * l[i];

  }
  phi_kakko = phi_kakko / L;

  V_kakko = 0.0;
  for( i = 1; i <= NUM; i++ ){

    V_kakko += kappa[i] * V[i];

  }
  V_kakko = phi_kakko / L;
  
  for( i = 1; i <= NUM; i++ ){

    psi[i] = phiphi(kappa[i]) * l[i] * ( ( f_kakko / phi_kakko ) - ( f[i] / phiphi(kappa[i]) ) + ( ( L * phi_kakko / ( NUM * l[i] * phiphi(kappa[i]) ) ) - 1 ) * omega(NUM) );
    //psi[i] = phiphi(kappa[i]) * l[i] * ( ( f_kakko / phi_kakko ) - ( f[i] / phiphi(kappa[i]) ) + ( ( L * phi_kakko / ( NUM * l[i] * phiphi(kappa[i]) ) ) - 1 ) * ( 0.1 / dt ) );

    //psi[i] = phiphi(kappa[i]) * l[i] * ( ( f_kakko / phi_kakko ) - ( f[i] / phiphi(kappa[i]) ) + ( ( L * phi_kakko / ( NUM * l[i] * phiphi(kappa[i]) ) ) - 1 ) * ( 100 + 100 * V_kakko ) );
    
  }

  for( i = 1; i <= NUM; i++ ){

    PSI[i] = 0.0;

  }

  for( i = 2; i <= NUM; i++ ){
    
    for( j = 2; j <= i; j++ ){
    
      PSI[i] += psi[j];
	
    }
    
  }
  

  W[1] = 0.0;
  for( i = 2; i <= NUM; i++ ){
    
    W[1] += PSI[i] * l_mid[i];
  
  }
  W[1] = -W[1] / ( L * phiphi(kappa_mid[1]) );
  
  for( i = 2; i <= NUM; i++ ){
    
    W[i] = ( phiphi(kappa_mid[1]) * W[1] + PSI[i] ) / phiphi(kappa_mid[i]);

  }
  connect(W);

  free(PSI);
  free(psi);
  free(f);

}

double phiphi( double x ){

  return( 1 - epsl + epsl * sqrt(1 - epsl + epsl * x * x) );
  
}

double phiphi_x( double x ){

  return( epsl * epsl * x / sqrt(1 - epsl + epsl * x * x) );
  
}

double Gamma( double l, int i, double *p1, double *p2 ){

  double theta;
  theta = theta0 + ( l * M_PI / 3.0 );

  return(
	 sqrt(
	      ( cos(theta) * p1[i] + sin(theta) * p2[i] )
	    * ( cos(theta) * p1[i] + sin(theta) * p2[i] )
	    + eps * eps
	    * ( -sin(theta) * p1[i] + cos(theta) * p2[i] )
	    * ( -sin(theta) * p1[i] + cos(theta) * p2[i] )
	     )
	);
  
}

double Gamma_1( double l, int i, double *p1, double *p2 ){

  double theta;
  theta = theta0 + ( l * M_PI / 3.0 );
  
  return(
	 cos(theta) * cos(theta) * p1[i] + cos(theta) * sin(theta) * p2[i]
	 + eps * eps * ( sin(theta) * sin(theta) * p1[i] - cos(theta) * sin(theta) * p2[i] )
	);
   
}

double Gamma_2( double l, int i, double *p1, double *p2 ){

  double theta;
  theta = theta0 + ( l * M_PI / 3.0 );
  
  return(
	 cos(theta) * sin(theta) * p1[i] + sin(theta) * sin(theta) * p2[i]
	 + eps * eps * ( -cos(theta) * sin(theta) * p1[i] + cos(theta) * cos(theta) * p2[i] )
	);
   
}

double Gamma_3( double l, int i ){

  double theta;
  theta = theta0 + ( l * M_PI / 3.0 );
  
  return( cos(theta) * cos(theta) + eps * eps * sin(theta) * sin(theta) );
  
}

double Gamma_4( double l, int i ){

  double theta;
  theta = theta0 + ( l * M_PI / 3.0 );
  
  return( sin(theta) * sin(theta) + eps * eps * cos(theta) * cos(theta) );

}

double Gamma_5( double l, int i ){

  double theta;
  theta = theta0 + ( l * M_PI / 3.0 );
  
  return( cos(theta) * sin(theta) - eps * eps * cos(theta) * sin(theta) );

}

double Gammap1( double l, int i, double *p1, double *p2 ){

  return( Gamma_1(l,i,p1,p2) / Gamma(l,i,p1,p2) );

}

double Gammap2( double l, int i, double *p1, double *p2 ){

  return( Gamma_2(l,i,p1,p2) / Gamma(l,i,p1,p2) );

}

double Gamma_p1p1( double l, int i, double *p1, double *p2 ){
  
  return( -( Gamma_1(l,i,p1,p2) * Gamma_1(l,i,p1,p2) / ( Gamma(l,i,p1,p2) * Gamma(l,i,p1,p2) * Gamma(l,i,p1,p2) ) ) + ( Gamma_3(l,i) / Gamma(l,i,p1,p2) ) );
  
}

double Gamma_p2p2( double l, int i, double *p1, double *p2 ){
  
  return( -( Gamma_2(l,i,p1,p2) * Gamma_2(l,i,p1,p2) / ( Gamma(l,i,p1,p2) * Gamma(l,i,p1,p2) * Gamma(l,i,p1,p2) ) ) + ( Gamma_4(l,i) / Gamma(l,i,p1,p2) ) );
  
}

double Gamma_p1p2( double l, int i, double *p1, double *p2 ){
  
  return( -( Gamma_1(l,i,p1,p2) * Gamma_2(l,i,p1,p2) / ( Gamma(l,i,p1,p2) * Gamma(l,i,p1,p2) * Gamma(l,i,p1,p2) ) ) + ( Gamma_5(l,i) / Gamma(l,i,p1,p2) ) );
  
}
