      SUBROUTINE CHET21( ITYPE, UPLO, N, KBAND, A, LDA, D, E, U, LDU, V, LDV, TAU, WORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                ITYPE, KBAND, LDA, LDU, LDV, N;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * ), RESULT( 2 ), RWORK( * )
      COMPLEX            A( LDA, * ), TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TEN
      const              ZERO = 0.0E+0, ONE = 1.0E+0, TEN = 10.0E+0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LOWER;
      String             CUPLO;
      int                IINFO, J, JCOL, JR, JROW;
      REAL               ANORM, ULP, UNFL, WNORM
      COMPLEX            VSAVE
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGE, CLANHE, SLAMCH
      // EXTERNAL LSAME, CLANGE, CLANHE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CHER, CHER2, CLACPY, CLARFY, CLASET, CUNM2L, CUNM2R
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      RESULT( 1 ) = ZERO
      if (ITYPE == 1) RESULT( 2 ) = ZERO       IF( N.LE.0 ) RETURN;

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
         ANORM = MAX( CLANHE( '1', CUPLO, N, A, LDA, RWORK ), UNFL )
      }

      // Compute error matrix:

      if ( ITYPE == 1 ) {

         // ITYPE=1: error = A - U S U**H

         claset('Full', N, N, CZERO, CZERO, WORK, N );
         clacpy(CUPLO, N, N, A, LDA, WORK, N );

         for (J = 1; J <= N; J++) { // 10
            cher(CUPLO, N, -D( J ), U( 1, J ), 1, WORK, N );
         } // 10

         if ( N.GT.1 && KBAND == 1 ) {
            for (J = 1; J <= N - 1; J++) { // 20
               cher2(CUPLO, N, -CMPLX( E( J ) ), U( 1, J ), 1, U( 1, J+1 ), 1, WORK, N );
            } // 20
         }
         WNORM = CLANHE( '1', CUPLO, N, WORK, N, RWORK )

      } else if ( ITYPE == 2 ) {

         // ITYPE=2: error = V S V**H - A

         claset('Full', N, N, CZERO, CZERO, WORK, N );

         if ( LOWER ) {
            WORK( N**2 ) = D( N )
            DO 40 J = N - 1, 1, -1
               if ( KBAND == 1 ) {
                  WORK( ( N+1 )*( J-1 )+2 ) = ( CONE-TAU( J ) )*E( J )
                  for (JR = J + 2; JR <= N; JR++) { // 30
                     WORK( ( J-1 )*N+JR ) = -TAU( J )*E( J )*V( JR, J )
                  } // 30
               }

               VSAVE = V( J+1, J )
               V( J+1, J ) = ONE
               clarfy('L', N-J, V( J+1, J ), 1, TAU( J ), WORK( ( N+1 )*J+1 ), N, WORK( N**2+1 ) );
               V( J+1, J ) = VSAVE
               WORK( ( N+1 )*( J-1 )+1 ) = D( J )
            } // 40
         } else {
            WORK( 1 ) = D( 1 )
            for (J = 1; J <= N - 1; J++) { // 60
               if ( KBAND == 1 ) {
                  WORK( ( N+1 )*J ) = ( CONE-TAU( J ) )*E( J )
                  for (JR = 1; JR <= J - 1; JR++) { // 50
                     WORK( J*N+JR ) = -TAU( J )*E( J )*V( JR, J+1 )
                  } // 50
               }

               VSAVE = V( J, J+1 )
               V( J, J+1 ) = ONE
               clarfy('U', J, V( 1, J+1 ), 1, TAU( J ), WORK, N, WORK( N**2+1 ) );
               V( J, J+1 ) = VSAVE
               WORK( ( N+1 )*J+1 ) = D( J+1 )
            } // 60
         }

         for (JCOL = 1; JCOL <= N; JCOL++) { // 90
            if ( LOWER ) {
               for (JROW = JCOL; JROW <= N; JROW++) { // 70
                  WORK( JROW+N*( JCOL-1 ) ) = WORK( JROW+N*( JCOL-1 ) ) - A( JROW, JCOL )
               } // 70
            } else {
               for (JROW = 1; JROW <= JCOL; JROW++) { // 80
                  WORK( JROW+N*( JCOL-1 ) ) = WORK( JROW+N*( JCOL-1 ) ) - A( JROW, JCOL )
               } // 80
            }
         } // 90
         WNORM = CLANHE( '1', CUPLO, N, WORK, N, RWORK )

      } else if ( ITYPE == 3 ) {

         // ITYPE=3: error = U V**H - I

         if (N.LT.2) RETURN;
         clacpy(' ', N, N, U, LDU, WORK, N );
         if ( LOWER ) {
            cunm2r('R', 'C', N, N-1, N-1, V( 2, 1 ), LDV, TAU, WORK( N+1 ), N, WORK( N**2+1 ), IINFO );
         } else {
            cunm2l('R', 'C', N, N-1, N-1, V( 1, 2 ), LDV, TAU, WORK, N, WORK( N**2+1 ), IINFO );
         }
         if ( IINFO != 0 ) {
            RESULT( 1 ) = TEN / ULP
            RETURN
         }

         for (J = 1; J <= N; J++) { // 100
            WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - CONE
         } // 100

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

         for (J = 1; J <= N; J++) { // 110
            WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - CONE
         } // 110

         RESULT( 2 ) = MIN( CLANGE( '1', N, N, WORK, N, RWORK ), REAL( N ) ) / ( N*ULP )
      }

      RETURN

      // End of CHET21

      }
