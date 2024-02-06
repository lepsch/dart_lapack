      void cspt03(UPLO, N, A, AINV, WORK, LDW, RWORK, RCOND, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                LDW, N;
      double               RCOND, RESID;
      double               RWORK( * );
      Complex            A( * ), AINV( * ), WORK( LDW, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, ICOL, J, JCOL, K, KCOL, NALL;
      double               AINVNM, ANORM, EPS;
      Complex            T;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANGE, CLANSP, SLAMCH;
      //- COMPLEX            CDOTU;
      // EXTERNAL lsame, CLANGE, CLANSP, SLAMCH, CDOTU
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RCOND = ONE;
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = CLANSP( '1', UPLO, N, A, RWORK );
      AINVNM = CLANSP( '1', UPLO, N, AINV, RWORK );
      if ( ANORM <= ZERO || AINVNM <= ZERO ) {
         RCOND = ZERO;
         RESID = ONE / EPS;
         return;
      }
      RCOND = ( ONE/ANORM ) / AINVNM;

      // Case where both A and AINV are upper triangular:
      // Each element of - A * AINV is computed by taking the dot product
      // of a row of A with a column of AINV.

      if ( lsame( UPLO, 'U' ) ) {
         for (I = 1; I <= N; I++) { // 70
            ICOL = ( ( I-1 )*I ) / 2 + 1;

            // Code when J <= I

            for (J = 1; J <= I; J++) { // 30
               JCOL = ( ( J-1 )*J ) / 2 + 1;
               T = CDOTU( J, A( ICOL ), 1, AINV( JCOL ), 1 );
               JCOL = JCOL + 2*J - 1;
               KCOL = ICOL - 1;
               for (K = J + 1; K <= I; K++) { // 10
                  T = T + A( KCOL+K )*AINV( JCOL );
                  JCOL = JCOL + K;
               } // 10
               KCOL = KCOL + 2*I;
               for (K = I + 1; K <= N; K++) { // 20
                  T = T + A( KCOL )*AINV( JCOL );
                  KCOL = KCOL + K;
                  JCOL = JCOL + K;
               } // 20
               WORK[I][J] = -T;
            } // 30

            // Code when J > I

            for (J = I + 1; J <= N; J++) { // 60
               JCOL = ( ( J-1 )*J ) / 2 + 1;
               T = CDOTU( I, A( ICOL ), 1, AINV( JCOL ), 1 );
               JCOL = JCOL - 1;
               KCOL = ICOL + 2*I - 1;
               for (K = I + 1; K <= J; K++) { // 40
                  T = T + A( KCOL )*AINV( JCOL+K );
                  KCOL = KCOL + K;
               } // 40
               JCOL = JCOL + 2*J;
               for (K = J + 1; K <= N; K++) { // 50
                  T = T + A( KCOL )*AINV( JCOL );
                  KCOL = KCOL + K;
                  JCOL = JCOL + K;
               } // 50
               WORK[I][J] = -T;
            } // 60
         } // 70
      } else {

         // Case where both A and AINV are lower triangular

         NALL = ( N*( N+1 ) ) / 2;
         for (I = 1; I <= N; I++) { // 140

            // Code when J <= I

            ICOL = NALL - ( ( N-I+1 )*( N-I+2 ) ) / 2 + 1;
            for (J = 1; J <= I; J++) { // 100
               JCOL = NALL - ( ( N-J )*( N-J+1 ) ) / 2 - ( N-I );
               T = CDOTU( N-I+1, A( ICOL ), 1, AINV( JCOL ), 1 );
               KCOL = I;
               JCOL = J;
               for (K = 1; K <= J - 1; K++) { // 80
                  T = T + A( KCOL )*AINV( JCOL );
                  JCOL = JCOL + N - K;
                  KCOL = KCOL + N - K;
               } // 80
               JCOL = JCOL - J;
               for (K = J; K <= I - 1; K++) { // 90
                  T = T + A( KCOL )*AINV( JCOL+K );
                  KCOL = KCOL + N - K;
               } // 90
               WORK[I][J] = -T;
            } // 100

            // Code when J > I

            ICOL = NALL - ( ( N-I )*( N-I+1 ) ) / 2;
            for (J = I + 1; J <= N; J++) { // 130
               JCOL = NALL - ( ( N-J+1 )*( N-J+2 ) ) / 2 + 1;
               T = CDOTU( N-J+1, A( ICOL-N+J ), 1, AINV( JCOL ), 1 );
               KCOL = I;
               JCOL = J;
               for (K = 1; K <= I - 1; K++) { // 110
                  T = T + A( KCOL )*AINV( JCOL );
                  JCOL = JCOL + N - K;
                  KCOL = KCOL + N - K;
               } // 110
               KCOL = KCOL - I;
               for (K = I; K <= J - 1; K++) { // 120
                  T = T + A( KCOL+K )*AINV( JCOL );
                  JCOL = JCOL + N - K;
               } // 120
               WORK[I][J] = -T;
            } // 130
         } // 140
      }

      // Add the identity matrix to WORK .

      for (I = 1; I <= N; I++) { // 150
         WORK[I][I] = WORK( I, I ) + ONE;
      } // 150

      // Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)

      RESID = CLANGE( '1', N, N, WORK, LDW, RWORK );

      RESID = ( ( RESID*RCOND )/EPS ) / REAL( N );

      return;
      }
