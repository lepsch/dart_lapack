      void slarfx(final int SIDE, final int M, final int N, final int V, final int TAU, final Matrix<double> C_, final int LDC, final Array<double> WORK_,) {
  final C = C_.dim();
  final WORK = WORK_.dim();

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             SIDE;
      int                LDC, M, N;
      double               TAU;
      double               C( LDC, * ), V( * ), WORK( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                J;
      double               SUM, T1, T10, T2, T3, T4, T5, T6, T7, T8, T9, V1, V10, V2, V3, V4, V5, V6, V7, V8, V9;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARF

      if (TAU == ZERO) return;
      if ( lsame( SIDE, 'L' ) ) {

         // Form  H * C, where H has order m.

         GO TO ( 10, 30, 50, 70, 90, 110, 130, 150, 170, 190 )M;

         // Code for general M

         slarf(SIDE, M, N, V, 1, TAU, C, LDC, WORK );
         GO TO 410;
         } // 10

         // Special code for 1 x 1 Householder

         T1 = ONE - TAU*V( 1 )*V( 1 );
         for (J = 1; J <= N; J++) { // 20
            C[1][J] = T1*C( 1, J );
         } // 20
         GO TO 410;
         } // 30

         // Special code for 2 x 2 Householder

         V1 = V( 1 );
         T1 = TAU*V1;
         V2 = V( 2 );
         T2 = TAU*V2;
         for (J = 1; J <= N; J++) { // 40
            SUM = V1*C( 1, J ) + V2*C( 2, J );
            C[1][J] = C( 1, J ) - SUM*T1;
            C[2][J] = C( 2, J ) - SUM*T2;
         } // 40
         GO TO 410;
         } // 50

         // Special code for 3 x 3 Householder

         V1 = V( 1 );
         T1 = TAU*V1;
         V2 = V( 2 );
         T2 = TAU*V2;
         V3 = V( 3 );
         T3 = TAU*V3;
         for (J = 1; J <= N; J++) { // 60
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J );
            C[1][J] = C( 1, J ) - SUM*T1;
            C[2][J] = C( 2, J ) - SUM*T2;
            C[3][J] = C( 3, J ) - SUM*T3;
         } // 60
         GO TO 410;
         } // 70

         // Special code for 4 x 4 Householder

         V1 = V( 1 );
         T1 = TAU*V1;
         V2 = V( 2 );
         T2 = TAU*V2;
         V3 = V( 3 );
         T3 = TAU*V3;
         V4 = V( 4 );
         T4 = TAU*V4;
         for (J = 1; J <= N; J++) { // 80
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + V4*C( 4, J );
            C[1][J] = C( 1, J ) - SUM*T1;
            C[2][J] = C( 2, J ) - SUM*T2;
            C[3][J] = C( 3, J ) - SUM*T3;
            C[4][J] = C( 4, J ) - SUM*T4;
         } // 80
         GO TO 410;
         } // 90

         // Special code for 5 x 5 Householder

         V1 = V( 1 );
         T1 = TAU*V1;
         V2 = V( 2 );
         T2 = TAU*V2;
         V3 = V( 3 );
         T3 = TAU*V3;
         V4 = V( 4 );
         T4 = TAU*V4;
         V5 = V( 5 );
         T5 = TAU*V5;
         for (J = 1; J <= N; J++) { // 100
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + V4*C( 4, J ) + V5*C( 5, J );
            C[1][J] = C( 1, J ) - SUM*T1;
            C[2][J] = C( 2, J ) - SUM*T2;
            C[3][J] = C( 3, J ) - SUM*T3;
            C[4][J] = C( 4, J ) - SUM*T4;
            C[5][J] = C( 5, J ) - SUM*T5;
         } // 100
         GO TO 410;
         } // 110

         // Special code for 6 x 6 Householder

         V1 = V( 1 );
         T1 = TAU*V1;
         V2 = V( 2 );
         T2 = TAU*V2;
         V3 = V( 3 );
         T3 = TAU*V3;
         V4 = V( 4 );
         T4 = TAU*V4;
         V5 = V( 5 );
         T5 = TAU*V5;
         V6 = V( 6 );
         T6 = TAU*V6;
         for (J = 1; J <= N; J++) { // 120
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J );
            C[1][J] = C( 1, J ) - SUM*T1;
            C[2][J] = C( 2, J ) - SUM*T2;
            C[3][J] = C( 3, J ) - SUM*T3;
            C[4][J] = C( 4, J ) - SUM*T4;
            C[5][J] = C( 5, J ) - SUM*T5;
            C[6][J] = C( 6, J ) - SUM*T6;
         } // 120
         GO TO 410;
         } // 130

         // Special code for 7 x 7 Householder

         V1 = V( 1 );
         T1 = TAU*V1;
         V2 = V( 2 );
         T2 = TAU*V2;
         V3 = V( 3 );
         T3 = TAU*V3;
         V4 = V( 4 );
         T4 = TAU*V4;
         V5 = V( 5 );
         T5 = TAU*V5;
         V6 = V( 6 );
         T6 = TAU*V6;
         V7 = V( 7 );
         T7 = TAU*V7;
         for (J = 1; J <= N; J++) { // 140
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + V7*C( 7, J );
            C[1][J] = C( 1, J ) - SUM*T1;
            C[2][J] = C( 2, J ) - SUM*T2;
            C[3][J] = C( 3, J ) - SUM*T3;
            C[4][J] = C( 4, J ) - SUM*T4;
            C[5][J] = C( 5, J ) - SUM*T5;
            C[6][J] = C( 6, J ) - SUM*T6;
            C[7][J] = C( 7, J ) - SUM*T7;
         } // 140
         GO TO 410;
         } // 150

         // Special code for 8 x 8 Householder

         V1 = V( 1 );
         T1 = TAU*V1;
         V2 = V( 2 );
         T2 = TAU*V2;
         V3 = V( 3 );
         T3 = TAU*V3;
         V4 = V( 4 );
         T4 = TAU*V4;
         V5 = V( 5 );
         T5 = TAU*V5;
         V6 = V( 6 );
         T6 = TAU*V6;
         V7 = V( 7 );
         T7 = TAU*V7;
         V8 = V( 8 );
         T8 = TAU*V8;
         for (J = 1; J <= N; J++) { // 160
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + V7*C( 7, J ) + V8*C( 8, J );
            C[1][J] = C( 1, J ) - SUM*T1;
            C[2][J] = C( 2, J ) - SUM*T2;
            C[3][J] = C( 3, J ) - SUM*T3;
            C[4][J] = C( 4, J ) - SUM*T4;
            C[5][J] = C( 5, J ) - SUM*T5;
            C[6][J] = C( 6, J ) - SUM*T6;
            C[7][J] = C( 7, J ) - SUM*T7;
            C[8][J] = C( 8, J ) - SUM*T8;
         } // 160
         GO TO 410;
         } // 170

         // Special code for 9 x 9 Householder

         V1 = V( 1 );
         T1 = TAU*V1;
         V2 = V( 2 );
         T2 = TAU*V2;
         V3 = V( 3 );
         T3 = TAU*V3;
         V4 = V( 4 );
         T4 = TAU*V4;
         V5 = V( 5 );
         T5 = TAU*V5;
         V6 = V( 6 );
         T6 = TAU*V6;
         V7 = V( 7 );
         T7 = TAU*V7;
         V8 = V( 8 );
         T8 = TAU*V8;
         V9 = V( 9 );
         T9 = TAU*V9;
         for (J = 1; J <= N; J++) { // 180
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + V7*C( 7, J ) + V8*C( 8, J ) + V9*C( 9, J );
            C[1][J] = C( 1, J ) - SUM*T1;
            C[2][J] = C( 2, J ) - SUM*T2;
            C[3][J] = C( 3, J ) - SUM*T3;
            C[4][J] = C( 4, J ) - SUM*T4;
            C[5][J] = C( 5, J ) - SUM*T5;
            C[6][J] = C( 6, J ) - SUM*T6;
            C[7][J] = C( 7, J ) - SUM*T7;
            C[8][J] = C( 8, J ) - SUM*T8;
            C[9][J] = C( 9, J ) - SUM*T9;
         } // 180
         GO TO 410;
         } // 190

         // Special code for 10 x 10 Householder

         V1 = V( 1 );
         T1 = TAU*V1;
         V2 = V( 2 );
         T2 = TAU*V2;
         V3 = V( 3 );
         T3 = TAU*V3;
         V4 = V( 4 );
         T4 = TAU*V4;
         V5 = V( 5 );
         T5 = TAU*V5;
         V6 = V( 6 );
         T6 = TAU*V6;
         V7 = V( 7 );
         T7 = TAU*V7;
         V8 = V( 8 );
         T8 = TAU*V8;
         V9 = V( 9 );
         T9 = TAU*V9;
         V10 = V( 10 );
         T10 = TAU*V10;
         for (J = 1; J <= N; J++) { // 200
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + V7*C( 7, J ) + V8*C( 8, J ) + V9*C( 9, J ) + V10*C( 10, J );
            C[1][J] = C( 1, J ) - SUM*T1;
            C[2][J] = C( 2, J ) - SUM*T2;
            C[3][J] = C( 3, J ) - SUM*T3;
            C[4][J] = C( 4, J ) - SUM*T4;
            C[5][J] = C( 5, J ) - SUM*T5;
            C[6][J] = C( 6, J ) - SUM*T6;
            C[7][J] = C( 7, J ) - SUM*T7;
            C[8][J] = C( 8, J ) - SUM*T8;
            C[9][J] = C( 9, J ) - SUM*T9;
            C[10][J] = C( 10, J ) - SUM*T10;
         } // 200
         GO TO 410;
      } else {

         // Form  C * H, where H has order n.

         GO TO ( 210, 230, 250, 270, 290, 310, 330, 350, 370, 390 )N;

         // Code for general N

         slarf(SIDE, M, N, V, 1, TAU, C, LDC, WORK );
         GO TO 410;
         } // 210

         // Special code for 1 x 1 Householder

         T1 = ONE - TAU*V( 1 )*V( 1 );
         for (J = 1; J <= M; J++) { // 220
            C[J][1] = T1*C( J, 1 );
         } // 220
         GO TO 410;
         } // 230

         // Special code for 2 x 2 Householder

         V1 = V( 1 );
         T1 = TAU*V1;
         V2 = V( 2 );
         T2 = TAU*V2;
         for (J = 1; J <= M; J++) { // 240
            SUM = V1*C( J, 1 ) + V2*C( J, 2 );
            C[J][1] = C( J, 1 ) - SUM*T1;
            C[J][2] = C( J, 2 ) - SUM*T2;
         } // 240
         GO TO 410;
         } // 250

         // Special code for 3 x 3 Householder

         V1 = V( 1 );
         T1 = TAU*V1;
         V2 = V( 2 );
         T2 = TAU*V2;
         V3 = V( 3 );
         T3 = TAU*V3;
         for (J = 1; J <= M; J++) { // 260
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 );
            C[J][1] = C( J, 1 ) - SUM*T1;
            C[J][2] = C( J, 2 ) - SUM*T2;
            C[J][3] = C( J, 3 ) - SUM*T3;
         } // 260
         GO TO 410;
         } // 270

         // Special code for 4 x 4 Householder

         V1 = V( 1 );
         T1 = TAU*V1;
         V2 = V( 2 );
         T2 = TAU*V2;
         V3 = V( 3 );
         T3 = TAU*V3;
         V4 = V( 4 );
         T4 = TAU*V4;
         for (J = 1; J <= M; J++) { // 280
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + V4*C( J, 4 );
            C[J][1] = C( J, 1 ) - SUM*T1;
            C[J][2] = C( J, 2 ) - SUM*T2;
            C[J][3] = C( J, 3 ) - SUM*T3;
            C[J][4] = C( J, 4 ) - SUM*T4;
         } // 280
         GO TO 410;
         } // 290

         // Special code for 5 x 5 Householder

         V1 = V( 1 );
         T1 = TAU*V1;
         V2 = V( 2 );
         T2 = TAU*V2;
         V3 = V( 3 );
         T3 = TAU*V3;
         V4 = V( 4 );
         T4 = TAU*V4;
         V5 = V( 5 );
         T5 = TAU*V5;
         for (J = 1; J <= M; J++) { // 300
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + V4*C( J, 4 ) + V5*C( J, 5 );
            C[J][1] = C( J, 1 ) - SUM*T1;
            C[J][2] = C( J, 2 ) - SUM*T2;
            C[J][3] = C( J, 3 ) - SUM*T3;
            C[J][4] = C( J, 4 ) - SUM*T4;
            C[J][5] = C( J, 5 ) - SUM*T5;
         } // 300
         GO TO 410;
         } // 310

         // Special code for 6 x 6 Householder

         V1 = V( 1 );
         T1 = TAU*V1;
         V2 = V( 2 );
         T2 = TAU*V2;
         V3 = V( 3 );
         T3 = TAU*V3;
         V4 = V( 4 );
         T4 = TAU*V4;
         V5 = V( 5 );
         T5 = TAU*V5;
         V6 = V( 6 );
         T6 = TAU*V6;
         for (J = 1; J <= M; J++) { // 320
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 );
            C[J][1] = C( J, 1 ) - SUM*T1;
            C[J][2] = C( J, 2 ) - SUM*T2;
            C[J][3] = C( J, 3 ) - SUM*T3;
            C[J][4] = C( J, 4 ) - SUM*T4;
            C[J][5] = C( J, 5 ) - SUM*T5;
            C[J][6] = C( J, 6 ) - SUM*T6;
         } // 320
         GO TO 410;
         } // 330

         // Special code for 7 x 7 Householder

         V1 = V( 1 );
         T1 = TAU*V1;
         V2 = V( 2 );
         T2 = TAU*V2;
         V3 = V( 3 );
         T3 = TAU*V3;
         V4 = V( 4 );
         T4 = TAU*V4;
         V5 = V( 5 );
         T5 = TAU*V5;
         V6 = V( 6 );
         T6 = TAU*V6;
         V7 = V( 7 );
         T7 = TAU*V7;
         for (J = 1; J <= M; J++) { // 340
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + V7*C( J, 7 );
            C[J][1] = C( J, 1 ) - SUM*T1;
            C[J][2] = C( J, 2 ) - SUM*T2;
            C[J][3] = C( J, 3 ) - SUM*T3;
            C[J][4] = C( J, 4 ) - SUM*T4;
            C[J][5] = C( J, 5 ) - SUM*T5;
            C[J][6] = C( J, 6 ) - SUM*T6;
            C[J][7] = C( J, 7 ) - SUM*T7;
         } // 340
         GO TO 410;
         } // 350

         // Special code for 8 x 8 Householder

         V1 = V( 1 );
         T1 = TAU*V1;
         V2 = V( 2 );
         T2 = TAU*V2;
         V3 = V( 3 );
         T3 = TAU*V3;
         V4 = V( 4 );
         T4 = TAU*V4;
         V5 = V( 5 );
         T5 = TAU*V5;
         V6 = V( 6 );
         T6 = TAU*V6;
         V7 = V( 7 );
         T7 = TAU*V7;
         V8 = V( 8 );
         T8 = TAU*V8;
         for (J = 1; J <= M; J++) { // 360
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + V7*C( J, 7 ) + V8*C( J, 8 );
            C[J][1] = C( J, 1 ) - SUM*T1;
            C[J][2] = C( J, 2 ) - SUM*T2;
            C[J][3] = C( J, 3 ) - SUM*T3;
            C[J][4] = C( J, 4 ) - SUM*T4;
            C[J][5] = C( J, 5 ) - SUM*T5;
            C[J][6] = C( J, 6 ) - SUM*T6;
            C[J][7] = C( J, 7 ) - SUM*T7;
            C[J][8] = C( J, 8 ) - SUM*T8;
         } // 360
         GO TO 410;
         } // 370

         // Special code for 9 x 9 Householder

         V1 = V( 1 );
         T1 = TAU*V1;
         V2 = V( 2 );
         T2 = TAU*V2;
         V3 = V( 3 );
         T3 = TAU*V3;
         V4 = V( 4 );
         T4 = TAU*V4;
         V5 = V( 5 );
         T5 = TAU*V5;
         V6 = V( 6 );
         T6 = TAU*V6;
         V7 = V( 7 );
         T7 = TAU*V7;
         V8 = V( 8 );
         T8 = TAU*V8;
         V9 = V( 9 );
         T9 = TAU*V9;
         for (J = 1; J <= M; J++) { // 380
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + V7*C( J, 7 ) + V8*C( J, 8 ) + V9*C( J, 9 );
            C[J][1] = C( J, 1 ) - SUM*T1;
            C[J][2] = C( J, 2 ) - SUM*T2;
            C[J][3] = C( J, 3 ) - SUM*T3;
            C[J][4] = C( J, 4 ) - SUM*T4;
            C[J][5] = C( J, 5 ) - SUM*T5;
            C[J][6] = C( J, 6 ) - SUM*T6;
            C[J][7] = C( J, 7 ) - SUM*T7;
            C[J][8] = C( J, 8 ) - SUM*T8;
            C[J][9] = C( J, 9 ) - SUM*T9;
         } // 380
         GO TO 410;
         } // 390

         // Special code for 10 x 10 Householder

         V1 = V( 1 );
         T1 = TAU*V1;
         V2 = V( 2 );
         T2 = TAU*V2;
         V3 = V( 3 );
         T3 = TAU*V3;
         V4 = V( 4 );
         T4 = TAU*V4;
         V5 = V( 5 );
         T5 = TAU*V5;
         V6 = V( 6 );
         T6 = TAU*V6;
         V7 = V( 7 );
         T7 = TAU*V7;
         V8 = V( 8 );
         T8 = TAU*V8;
         V9 = V( 9 );
         T9 = TAU*V9;
         V10 = V( 10 );
         T10 = TAU*V10;
         for (J = 1; J <= M; J++) { // 400
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + V7*C( J, 7 ) + V8*C( J, 8 ) + V9*C( J, 9 ) + V10*C( J, 10 );
            C[J][1] = C( J, 1 ) - SUM*T1;
            C[J][2] = C( J, 2 ) - SUM*T2;
            C[J][3] = C( J, 3 ) - SUM*T3;
            C[J][4] = C( J, 4 ) - SUM*T4;
            C[J][5] = C( J, 5 ) - SUM*T5;
            C[J][6] = C( J, 6 ) - SUM*T6;
            C[J][7] = C( J, 7 ) - SUM*T7;
            C[J][8] = C( J, 8 ) - SUM*T8;
            C[J][9] = C( J, 9 ) - SUM*T9;
            C[J][10] = C( J, 10 ) - SUM*T10;
         } // 400
         GO TO 410;
      }
  410 }
