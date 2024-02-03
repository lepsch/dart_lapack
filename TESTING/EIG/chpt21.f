      SUBROUTINE CHPT21( ITYPE, UPLO, N, KBAND, AP, D, E, U, LDU, VP, TAU, WORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                ITYPE, KBAND, LDU, N;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * ), RESULT( 2 ), RWORK( * )
      COMPLEX            AP( * ), TAU( * ), U( LDU, * ), VP( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TEN
      const              ZERO = 0.0E+0, ONE = 1.0E+0, TEN = 10.0E+0 ;
      REAL               HALF
      const              HALF = 1.0E+0 / 2.0E+0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LOWER;
      String             CUPLO;
      int                IINFO, J, JP, JP1, JR, LAP;
      REAL               ANORM, ULP, UNFL, WNORM
      COMPLEX            TEMP, VSAVE
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGE, CLANHP, SLAMCH
      COMPLEX            CDOTC
      // EXTERNAL LSAME, CLANGE, CLANHP, SLAMCH, CDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CGEMM, CHPMV, CHPR, CHPR2, CLACPY, CLASET, CUPMTR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // Constants

      RESULT( 1 ) = ZERO
      if (ITYPE == 1) RESULT( 2 ) = ZERO       IF( N.LE.0 ) RETURN;

      LAP = ( N*( N+1 ) ) / 2

      if ( LSAME( UPLO, 'U' ) ) {
         LOWER = false;
         CUPLO = 'U'
      } else {
         LOWER = true;
         CUPLO = 'L'
      }

      UNFL = SLAMCH( 'Safe minimum' )
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )

      // Some Error Checks

      if ( ITYPE.LT.1 || ITYPE.GT.3 ) {
         RESULT( 1 ) = TEN / ULP
         RETURN
      }

      // Do Test 1

      // Norm of A:

      if ( ITYPE == 3 ) {
         ANORM = ONE
      } else {
         ANORM = MAX( CLANHP( '1', CUPLO, N, AP, RWORK ), UNFL )
      }

      // Compute error matrix:

      if ( ITYPE == 1 ) {

         // ITYPE=1: error = A - U S U**H

         claset('Full', N, N, CZERO, CZERO, WORK, N );
         ccopy(LAP, AP, 1, WORK, 1 );

         for (J = 1; J <= N; J++) { // 10
            chpr(CUPLO, N, -D( J ), U( 1, J ), 1, WORK );
         } // 10

         if ( N.GT.1 && KBAND == 1 ) {
            for (J = 2; J <= N - 1; J++) { // 20
               chpr2(CUPLO, N, -CMPLX( E( J ) ), U( 1, J ), 1, U( 1, J-1 ), 1, WORK );
            } // 20
         }
         WNORM = CLANHP( '1', CUPLO, N, WORK, RWORK )

      } else if ( ITYPE == 2 ) {

         // ITYPE=2: error = V S V**H - A

         claset('Full', N, N, CZERO, CZERO, WORK, N );

         if ( LOWER ) {
            WORK( LAP ) = D( N )
            DO 40 J = N - 1, 1, -1
               JP = ( ( 2*N-J )*( J-1 ) ) / 2
               JP1 = JP + N - J
               if ( KBAND == 1 ) {
                  WORK( JP+J+1 ) = ( CONE-TAU( J ) )*E( J )
                  for (JR = J + 2; JR <= N; JR++) { // 30
                     WORK( JP+JR ) = -TAU( J )*E( J )*VP( JP+JR )
                  } // 30
               }

               if ( TAU( J ) != CZERO ) {
                  VSAVE = VP( JP+J+1 )
                  VP( JP+J+1 ) = CONE
                  chpmv('L', N-J, CONE, WORK( JP1+J+1 ), VP( JP+J+1 ), 1, CZERO, WORK( LAP+1 ), 1 )                   TEMP = -HALF*TAU( J )*CDOTC( N-J, WORK( LAP+1 ), 1, VP( JP+J+1 ), 1 );
                  caxpy(N-J, TEMP, VP( JP+J+1 ), 1, WORK( LAP+1 ), 1 );
                  chpr2('L', N-J, -TAU( J ), VP( JP+J+1 ), 1, WORK( LAP+1 ), 1, WORK( JP1+J+1 ) );

                  VP( JP+J+1 ) = VSAVE
               }
               WORK( JP+J ) = D( J )
            } // 40
         } else {
            WORK( 1 ) = D( 1 )
            for (J = 1; J <= N - 1; J++) { // 60
               JP = ( J*( J-1 ) ) / 2
               JP1 = JP + J
               if ( KBAND == 1 ) {
                  WORK( JP1+J ) = ( CONE-TAU( J ) )*E( J )
                  for (JR = 1; JR <= J - 1; JR++) { // 50
                     WORK( JP1+JR ) = -TAU( J )*E( J )*VP( JP1+JR )
                  } // 50
               }

               if ( TAU( J ) != CZERO ) {
                  VSAVE = VP( JP1+J )
                  VP( JP1+J ) = CONE
                  chpmv('U', J, CONE, WORK, VP( JP1+1 ), 1, CZERO, WORK( LAP+1 ), 1 )                   TEMP = -HALF*TAU( J )*CDOTC( J, WORK( LAP+1 ), 1, VP( JP1+1 ), 1 );
                  caxpy(J, TEMP, VP( JP1+1 ), 1, WORK( LAP+1 ), 1 );
                  chpr2('U', J, -TAU( J ), VP( JP1+1 ), 1, WORK( LAP+1 ), 1, WORK );
                  VP( JP1+J ) = VSAVE
               }
               WORK( JP1+J+1 ) = D( J+1 )
            } // 60
         }

         for (J = 1; J <= LAP; J++) { // 70
            WORK( J ) = WORK( J ) - AP( J )
         } // 70
         WNORM = CLANHP( '1', CUPLO, N, WORK, RWORK )

      } else if ( ITYPE == 3 ) {

         // ITYPE=3: error = U V**H - I

         if (N.LT.2) RETURN;
         clacpy(' ', N, N, U, LDU, WORK, N );
         cupmtr('R', CUPLO, 'C', N, N, VP, TAU, WORK, N, WORK( N**2+1 ), IINFO );
         if ( IINFO != 0 ) {
            RESULT( 1 ) = TEN / ULP
            RETURN
         }

         for (J = 1; J <= N; J++) { // 80
            WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - CONE
         } // 80

         WNORM = CLANGE( '1', N, N, WORK, N, RWORK )
      }

      if ( ANORM.GT.WNORM ) {
         RESULT( 1 ) = ( WNORM / ANORM ) / ( N*ULP )
      } else {
         if ( ANORM.LT.ONE ) {
            RESULT( 1 ) = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP )
         } else {
            RESULT( 1 ) = MIN( WNORM / ANORM, REAL( N ) ) / ( N*ULP )
         }
      }

      // Do Test 2

      // Compute  U U**H - I

      if ( ITYPE == 1 ) {
         cgemm('N', 'C', N, N, N, CONE, U, LDU, U, LDU, CZERO, WORK, N );

         for (J = 1; J <= N; J++) { // 90
            WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - CONE
         } // 90

         RESULT( 2 ) = MIN( CLANGE( '1', N, N, WORK, N, RWORK ), REAL( N ) ) / ( N*ULP )
      }

      RETURN

      // End of CHPT21

      }
