      void dspt21(ITYPE, UPLO, N, KBAND, AP, D, E, U, LDU, VP, TAU, WORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                ITYPE, KBAND, LDU, N;
      // ..
      // .. Array Arguments ..
      double             AP( * ), D( * ), E( * ), RESULT( 2 ), TAU( * ), U( LDU, * ), VP( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TEN;
      const              ZERO = 0.0, ONE = 1.0, TEN = 10.0 ;
      double             HALF;
      const              HALF = 1.0 / 2.0 ;
      // ..
      // .. Local Scalars ..
      bool               LOWER;
      String             CUPLO;
      int                IINFO, J, JP, JP1, JR, LAP;
      double             ANORM, TEMP, ULP, UNFL, VSAVE, WNORM;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- double             DDOT, DLAMCH, DLANGE, DLANSP;
      // EXTERNAL LSAME, DDOT, DLAMCH, DLANGE, DLANSP
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DGEMM, DLACPY, DLASET, DOPMTR, DSPMV, DSPR, DSPR2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // 1)      Constants

      RESULT( 1 ) = ZERO;
      if (ITYPE == 1) RESULT( 2 ) = ZERO;
      IF( N <= 0 ) return;

      LAP = ( N*( N+1 ) ) / 2;

      if ( LSAME( UPLO, 'U' ) ) {
         LOWER = false;
         CUPLO = 'U';
      } else {
         LOWER = true;
         CUPLO = 'L';
      }

      UNFL = DLAMCH( 'Safe minimum' );
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' );

      // Some Error Checks

      if ( ITYPE < 1 || ITYPE > 3 ) {
         RESULT( 1 ) = TEN / ULP;
         return;
      }

      // Do Test 1

      // Norm of A:

      if ( ITYPE == 3 ) {
         ANORM = ONE;
      } else {
         ANORM = max( DLANSP( '1', CUPLO, N, AP, WORK ), UNFL );
      }

      // Compute error matrix:

      if ( ITYPE == 1 ) {

         // ITYPE=1: error = A - U S U**T

         dlaset('Full', N, N, ZERO, ZERO, WORK, N );
         dcopy(LAP, AP, 1, WORK, 1 );

         for (J = 1; J <= N; J++) { // 10
            dspr(CUPLO, N, -D( J ), U( 1, J ), 1, WORK );
         } // 10

         if ( N > 1 && KBAND == 1 ) {
            for (J = 1; J <= N - 1; J++) { // 20
               dspr2(CUPLO, N, -E( J ), U( 1, J ), 1, U( 1, J+1 ), 1, WORK );
            } // 20
         }
         WNORM = DLANSP( '1', CUPLO, N, WORK, WORK( N**2+1 ) );

      } else if ( ITYPE == 2 ) {

         // ITYPE=2: error = V S V**T - A

         dlaset('Full', N, N, ZERO, ZERO, WORK, N );

         if ( LOWER ) {
            WORK( LAP ) = D( N );
            DO 40 J = N - 1, 1, -1;
               JP = ( ( 2*N-J )*( J-1 ) ) / 2;
               JP1 = JP + N - J;
               if ( KBAND == 1 ) {
                  WORK( JP+J+1 ) = ( ONE-TAU( J ) )*E( J );
                  for (JR = J + 2; JR <= N; JR++) { // 30
                     WORK( JP+JR ) = -TAU( J )*E( J )*VP( JP+JR );
                  } // 30
               }

               if ( TAU( J ) != ZERO ) {
                  VSAVE = VP( JP+J+1 );
                  VP( JP+J+1 ) = ONE;
                  dspmv('L', N-J, ONE, WORK( JP1+J+1 ), VP( JP+J+1 ), 1, ZERO, WORK( LAP+1 ), 1 )                   TEMP = -HALF*TAU( J )*DDOT( N-J, WORK( LAP+1 ), 1, VP( JP+J+1 ), 1 );
                  daxpy(N-J, TEMP, VP( JP+J+1 ), 1, WORK( LAP+1 ), 1 );
                  dspr2('L', N-J, -TAU( J ), VP( JP+J+1 ), 1, WORK( LAP+1 ), 1, WORK( JP1+J+1 ) );
                  VP( JP+J+1 ) = VSAVE;
               }
               WORK( JP+J ) = D( J );
            } // 40
         } else {
            WORK( 1 ) = D( 1 );
            for (J = 1; J <= N - 1; J++) { // 60
               JP = ( J*( J-1 ) ) / 2;
               JP1 = JP + J;
               if ( KBAND == 1 ) {
                  WORK( JP1+J ) = ( ONE-TAU( J ) )*E( J );
                  for (JR = 1; JR <= J - 1; JR++) { // 50
                     WORK( JP1+JR ) = -TAU( J )*E( J )*VP( JP1+JR );
                  } // 50
               }

               if ( TAU( J ) != ZERO ) {
                  VSAVE = VP( JP1+J );
                  VP( JP1+J ) = ONE;
                  dspmv('U', J, ONE, WORK, VP( JP1+1 ), 1, ZERO, WORK( LAP+1 ), 1 )                   TEMP = -HALF*TAU( J )*DDOT( J, WORK( LAP+1 ), 1, VP( JP1+1 ), 1 );
                  daxpy(J, TEMP, VP( JP1+1 ), 1, WORK( LAP+1 ), 1 );
                  dspr2('U', J, -TAU( J ), VP( JP1+1 ), 1, WORK( LAP+1 ), 1, WORK );
                  VP( JP1+J ) = VSAVE;
               }
               WORK( JP1+J+1 ) = D( J+1 );
            } // 60
         }

         for (J = 1; J <= LAP; J++) { // 70
            WORK( J ) = WORK( J ) - AP( J );
         } // 70
         WNORM = DLANSP( '1', CUPLO, N, WORK, WORK( LAP+1 ) );

      } else if ( ITYPE == 3 ) {

         // ITYPE=3: error = U V**T - I

         if (N < 2) return;
         dlacpy(' ', N, N, U, LDU, WORK, N );
         dopmtr('R', CUPLO, 'T', N, N, VP, TAU, WORK, N, WORK( N**2+1 ), IINFO );
         if ( IINFO != 0 ) {
            RESULT( 1 ) = TEN / ULP;
            return;
         }

         for (J = 1; J <= N; J++) { // 80
            WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - ONE;
         } // 80

         WNORM = DLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) );
      }

      if ( ANORM > WNORM ) {
         RESULT( 1 ) = ( WNORM / ANORM ) / ( N*ULP );
      } else {
         if ( ANORM < ONE ) {
            RESULT( 1 ) = ( min( WNORM, N*ANORM ) / ANORM ) / ( N*ULP );
         } else {
            RESULT( 1 ) = min( WNORM / ANORM, DBLE( N ) ) / ( N*ULP );
         }
      }

      // Do Test 2

      // Compute  U U**T - I

      if ( ITYPE == 1 ) {
         dgemm('N', 'C', N, N, N, ONE, U, LDU, U, LDU, ZERO, WORK, N );

         for (J = 1; J <= N; J++) { // 90
            WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - ONE;
         } // 90

         RESULT( 2 ) = min( DLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) ), DBLE( N ) ) / ( N*ULP );
      }

      return;
      }
