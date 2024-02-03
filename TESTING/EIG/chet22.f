      SUBROUTINE CHET22( ITYPE, UPLO, N, M, KBAND, A, LDA, D, E, U, LDU, V, LDV, TAU, WORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                ITYPE, KBAND, LDA, LDU, LDV, M, N;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * ), RESULT( 2 ), RWORK( * )
      COMPLEX            A( LDA, * ), TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E0, 0.0E0 ), CONE = ( 1.0E0, 0.0E0 ) ;
      // ..
      // .. Local Scalars ..
      int                J, JJ, JJ1, JJ2, NN, NNP1;
      REAL               ANORM, ULP, UNFL, WNORM
      // ..
      // .. External Functions ..
      REAL               CLANHE, SLAMCH
      // EXTERNAL CLANHE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CHEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      if (N.LE.0 .OR. M.LE.0) RETURN;

      UNFL = SLAMCH( 'Safe minimum' )
      ULP = SLAMCH( 'Precision' )

      // Do Test 1

      // Norm of A:

      ANORM = MAX( CLANHE( '1', UPLO, N, A, LDA, RWORK ), UNFL )

      // Compute error matrix:

      // ITYPE=1: error = U**H A U - S

      chemm('L', UPLO, N, M, CONE, A, LDA, U, LDU, CZERO, WORK, N );
      NN = N*N
      NNP1 = NN + 1
      cgemm('C', 'N', M, M, N, CONE, U, LDU, WORK, N, CZERO, WORK( NNP1 ), N );
      for (J = 1; J <= M; J++) { // 10
         JJ = NN + ( J-1 )*N + J
         WORK( JJ ) = WORK( JJ ) - D( J )
      } // 10
      if ( KBAND.EQ.1 .AND. N.GT.1 ) {
         for (J = 2; J <= M; J++) { // 20
            JJ1 = NN + ( J-1 )*N + J - 1
            JJ2 = NN + ( J-2 )*N + J
            WORK( JJ1 ) = WORK( JJ1 ) - E( J-1 )
            WORK( JJ2 ) = WORK( JJ2 ) - E( J-1 )
         } // 20
      }
      WNORM = CLANHE( '1', UPLO, M, WORK( NNP1 ), N, RWORK )

      if ( ANORM.GT.WNORM ) {
         RESULT( 1 ) = ( WNORM / ANORM ) / ( M*ULP )
      } else {
         if ( ANORM.LT.ONE ) {
            RESULT( 1 ) = ( MIN( WNORM, M*ANORM ) / ANORM ) / ( M*ULP )
         } else {
            RESULT( 1 ) = MIN( WNORM / ANORM, REAL( M ) ) / ( M*ULP )
         }
      }

      // Do Test 2

      // Compute  U**H U - I

      if (ITYPE.EQ.1) CALL CUNT01( 'Columns', N, M, U, LDU, WORK, 2*N*N, RWORK, RESULT( 2 ) );

      RETURN

      // End of CHET22

      }
