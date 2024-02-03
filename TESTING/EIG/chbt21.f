      SUBROUTINE CHBT21( UPLO, N, KA, KS, A, LDA, D, E, U, LDU, WORK, RWORK, RESULT );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                KA, KS, LDA, LDU, N;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * ), RESULT( 2 ), RWORK( * );
      COMPLEX            A( LDA, * ), U( LDU, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               LOWER;
      String             CUPLO;
      int                IKA, J, JC, JR;
      REAL               ANORM, ULP, UNFL, WNORM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGE, CLANHB, CLANHP, SLAMCH;
      // EXTERNAL LSAME, CLANGE, CLANHB, CLANHP, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CHPR, CHPR2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // Constants

      RESULT( 1 ) = ZERO;
      RESULT( 2 ) = ZERO;
      if (N <= 0) RETURN;

      IKA = MAX( 0, MIN( N-1, KA ) );

      if ( LSAME( UPLO, 'U' ) ) {
         LOWER = false;
         CUPLO = 'U';
      } else {
         LOWER = true;
         CUPLO = 'L';
      }

      UNFL = SLAMCH( 'Safe minimum' );
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' );

      // Some Error Checks

      // Do Test 1

      // Norm of A:

      ANORM = MAX( CLANHB( '1', CUPLO, N, IKA, A, LDA, RWORK ), UNFL );

      // Compute error matrix:    Error = A - U S U**H

      // Copy A from SB to SP storage format.

      J = 0;
      for (JC = 1; JC <= N; JC++) { // 50
         if ( LOWER ) {
            DO 10 JR = 1, MIN( IKA+1, N+1-JC );
               J = J + 1;
               WORK( J ) = A( JR, JC );
            } // 10
            for (JR = IKA + 2; JR <= N + 1 - JC; JR++) { // 20
               J = J + 1;
               WORK( J ) = ZERO;
            } // 20
         } else {
            for (JR = IKA + 2; JR <= JC; JR++) { // 30
               J = J + 1;
               WORK( J ) = ZERO;
            } // 30
            DO 40 JR = MIN( IKA, JC-1 ), 0, -1;
               J = J + 1;
               WORK( J ) = A( IKA+1-JR, JC );
            } // 40
         }
      } // 50

      for (J = 1; J <= N; J++) { // 60
         chpr(CUPLO, N, -D( J ), U( 1, J ), 1, WORK );
      } // 60

      if ( N > 1 && KS == 1 ) {
         for (J = 1; J <= N - 1; J++) { // 70
            chpr2(CUPLO, N, -CMPLX( E( J ) ), U( 1, J ), 1, U( 1, J+1 ), 1, WORK );
         } // 70
      }
      WNORM = CLANHP( '1', CUPLO, N, WORK, RWORK );

      if ( ANORM > WNORM ) {
         RESULT( 1 ) = ( WNORM / ANORM ) / ( N*ULP );
      } else {
         if ( ANORM < ONE ) {
            RESULT( 1 ) = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP );
         } else {
            RESULT( 1 ) = MIN( WNORM / ANORM, REAL( N ) ) / ( N*ULP );
         }
      }

      // Do Test 2

      // Compute  U U**H - I

      cgemm('N', 'C', N, N, N, CONE, U, LDU, U, LDU, CZERO, WORK, N );

      for (J = 1; J <= N; J++) { // 80
         WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - CONE;
      } // 80

      RESULT( 2 ) = MIN( CLANGE( '1', N, N, WORK, N, RWORK ), REAL( N ) ) / ( N*ULP );

      return;

      // End of CHBT21

      }
