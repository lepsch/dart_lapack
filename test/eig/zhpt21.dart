      void zhpt21(final int ITYPE, final int UPLO, final int N, final int KBAND, final int AP, final int D, final int E, final Matrix<double> U, final int LDU, final int VP, final int TAU, final Array<double> _WORK, final Array<double> RWORK, final int RESULT,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                ITYPE, KBAND, LDU, N;
      double             D( * ), E( * ), RESULT( 2 ), RWORK( * );
      Complex         AP( * ), TAU( * ), U( LDU, * ), VP( * ), WORK( * );
      // ..

      double             ZERO, ONE, TEN;
      const              ZERO = 0.0, ONE = 1.0, TEN = 10.0 ;
      double             HALF;
      const              HALF = 1.0 / 2.0 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      bool               LOWER;
      String             CUPLO;
      int                IINFO, J, JP, JP1, JR, LAP;
      double             ANORM, ULP, UNFL, WNORM;
      Complex         TEMP, VSAVE;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, ZLANGE, ZLANHP;
      //- Complex         ZDOTC;
      // EXTERNAL lsame, DLAMCH, ZLANGE, ZLANHP, ZDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZAXPY, ZCOPY, ZGEMM, ZHPMV, ZHPR, ZHPR2, ZLACPY, ZLASET, ZUPMTR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX, MIN

      // Constants

      RESULT[1] = ZERO;
      if (ITYPE == 1) RESULT( 2 ) = ZERO;
      IF( N <= 0 ) return;

      LAP = ( N*( N+1 ) ) / 2;

      if ( lsame( UPLO, 'U' ) ) {
         LOWER = false;
         CUPLO = 'U';
      } else {
         LOWER = true;
         CUPLO = 'L';
      }

      UNFL = dlamch( 'Safe minimum' );
      ULP = dlamch( 'Epsilon' )*dlamch( 'Base' );

      // Some Error Checks

      if ( ITYPE < 1 || ITYPE > 3 ) {
         RESULT[1] = TEN / ULP;
         return;
      }

      // Do Test 1

      // Norm of A:

      if ( ITYPE == 3 ) {
         ANORM = ONE;
      } else {
         ANORM = max( ZLANHP( '1', CUPLO, N, AP, RWORK ), UNFL );
      }

      // Compute error matrix:

      if ( ITYPE == 1 ) {

         // ITYPE=1: error = A - U S U**H

         zlaset('Full', N, N, CZERO, CZERO, WORK, N );
         zcopy(LAP, AP, 1, WORK, 1 );

         for (J = 1; J <= N; J++) { // 10
            zhpr(CUPLO, N, -D( J ), U( 1, J ), 1, WORK );
         } // 10

         if ( N > 1 && KBAND == 1 ) {
            for (J = 2; J <= N - 1; J++) { // 20
               zhpr2(CUPLO, N, -DCMPLX( E( J ) ), U( 1, J ), 1, U( 1, J-1 ), 1, WORK );
            } // 20
         }
         WNORM = ZLANHP( '1', CUPLO, N, WORK, RWORK );

      } else if ( ITYPE == 2 ) {

         // ITYPE=2: error = V S V**H - A

         zlaset('Full', N, N, CZERO, CZERO, WORK, N );

         if ( LOWER ) {
            WORK[LAP] = D( N );
            for (J = N - 1; J >= 1; J--) { // 40
               JP = ( ( 2*N-J )*( J-1 ) ) / 2;
               JP1 = JP + N - J;
               if ( KBAND == 1 ) {
                  WORK[JP+J+1] = ( CONE-TAU( J ) )*E( J );
                  for (JR = J + 2; JR <= N; JR++) { // 30
                     WORK[JP+JR] = -TAU( J )*E( J )*VP( JP+JR );
                  } // 30
               }

               if ( TAU( J ) != CZERO ) {
                  VSAVE = VP( JP+J+1 );
                  VP[JP+J+1] = CONE;
                  zhpmv('L', N-J, CONE, WORK( JP1+J+1 ), VP( JP+J+1 ), 1, CZERO, WORK( LAP+1 ), 1 )                   TEMP = -HALF*TAU( J )*ZDOTC( N-J, WORK( LAP+1 ), 1, VP( JP+J+1 ), 1 );
                  zaxpy(N-J, TEMP, VP( JP+J+1 ), 1, WORK( LAP+1 ), 1 );
                  zhpr2('L', N-J, -TAU( J ), VP( JP+J+1 ), 1, WORK( LAP+1 ), 1, WORK( JP1+J+1 ) );

                  VP[JP+J+1] = VSAVE;
               }
               WORK[JP+J] = D( J );
            } // 40
         } else {
            WORK[1] = D( 1 );
            for (J = 1; J <= N - 1; J++) { // 60
               JP = ( J*( J-1 ) ) / 2;
               JP1 = JP + J;
               if ( KBAND == 1 ) {
                  WORK[JP1+J] = ( CONE-TAU( J ) )*E( J );
                  for (JR = 1; JR <= J - 1; JR++) { // 50
                     WORK[JP1+JR] = -TAU( J )*E( J )*VP( JP1+JR );
                  } // 50
               }

               if ( TAU( J ) != CZERO ) {
                  VSAVE = VP( JP1+J );
                  VP[JP1+J] = CONE;
                  zhpmv('U', J, CONE, WORK, VP( JP1+1 ), 1, CZERO, WORK( LAP+1 ), 1 )                   TEMP = -HALF*TAU( J )*ZDOTC( J, WORK( LAP+1 ), 1, VP( JP1+1 ), 1 );
                  zaxpy(J, TEMP, VP( JP1+1 ), 1, WORK( LAP+1 ), 1 );
                  zhpr2('U', J, -TAU( J ), VP( JP1+1 ), 1, WORK( LAP+1 ), 1, WORK );
                  VP[JP1+J] = VSAVE;
               }
               WORK[JP1+J+1] = D( J+1 );
            } // 60
         }

         for (J = 1; J <= LAP; J++) { // 70
            WORK[J] = WORK( J ) - AP( J );
         } // 70
         WNORM = ZLANHP( '1', CUPLO, N, WORK, RWORK );

      } else if ( ITYPE == 3 ) {

         // ITYPE=3: error = U V**H - I

         if (N < 2) return;
         zlacpy(' ', N, N, U, LDU, WORK, N );
         zupmtr('R', CUPLO, 'C', N, N, VP, TAU, WORK, N, WORK( N**2+1 ), IINFO );
         if ( IINFO != 0 ) {
            RESULT[1] = TEN / ULP;
            return;
         }

         for (J = 1; J <= N; J++) { // 80
            WORK[( N+1 )*( J-1 )+1] = WORK( ( N+1 )*( J-1 )+1 ) - CONE;
         } // 80

         WNORM = ZLANGE( '1', N, N, WORK, N, RWORK );
      }

      if ( ANORM > WNORM ) {
         RESULT[1] = ( WNORM / ANORM ) / ( N*ULP );
      } else {
         if ( ANORM < ONE ) {
            RESULT[1] = ( min( WNORM, N*ANORM ) / ANORM ) / ( N*ULP );
         } else {
            RESULT[1] = min( WNORM / ANORM, N.toDouble() ) / ( N*ULP );
         }
      }

      // Do Test 2

      // Compute  U U**H - I

      if ( ITYPE == 1 ) {
         zgemm('N', 'C', N, N, N, CONE, U, LDU, U, LDU, CZERO, WORK, N );

         for (J = 1; J <= N; J++) { // 90
            WORK[( N+1 )*( J-1 )+1] = WORK( ( N+1 )*( J-1 )+1 ) - CONE;
         } // 90

         RESULT[2] = min( ZLANGE( '1', N, N, WORK, N, RWORK ), N.toDouble() ) / ( N*ULP );
      }

      }
