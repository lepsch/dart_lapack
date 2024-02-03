      SUBROUTINE DSBT21( UPLO, N, KA, KS, A, LDA, D, E, U, LDU, WORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                KA, KS, LDA, LDU, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), D( * ), E( * ), RESULT( 2 ), U( LDU, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      bool               LOWER;
      String             CUPLO;
      int                IKA, J, JC, JR, LW;
      double             ANORM, ULP, UNFL, WNORM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANGE, DLANSB, DLANSP;
      // EXTERNAL LSAME, DLAMCH, DLANGE, DLANSB, DLANSP
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DSPR, DSPR2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Constants

      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      if (N.LE.0) RETURN;

      IKA = MAX( 0, MIN( N-1, KA ) )
      LW = ( N*( N+1 ) ) / 2

      if ( LSAME( UPLO, 'U' ) ) {
         LOWER = false;
         CUPLO = 'U'
      } else {
         LOWER = true;
         CUPLO = 'L'
      }

      UNFL = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )

      // Some Error Checks

      // Do Test 1

      // Norm of A:

      ANORM = MAX( DLANSB( '1', CUPLO, N, IKA, A, LDA, WORK ), UNFL )

      // Compute error matrix:    Error = A - U S U**T

      // Copy A from SB to SP storage format.

      J = 0
      for (JC = 1; JC <= N; JC++) { // 50
         if ( LOWER ) {
            DO 10 JR = 1, MIN( IKA+1, N+1-JC )
               J = J + 1
               WORK( J ) = A( JR, JC )
            } // 10
            for (JR = IKA + 2; JR <= N + 1 - JC; JR++) { // 20
               J = J + 1
               WORK( J ) = ZERO
            } // 20
         } else {
            for (JR = IKA + 2; JR <= JC; JR++) { // 30
               J = J + 1
               WORK( J ) = ZERO
            } // 30
            DO 40 JR = MIN( IKA, JC-1 ), 0, -1
               J = J + 1
               WORK( J ) = A( IKA+1-JR, JC )
            } // 40
         }
      } // 50

      for (J = 1; J <= N; J++) { // 60
         dspr(CUPLO, N, -D( J ), U( 1, J ), 1, WORK );
      } // 60

      if ( N.GT.1 .AND. KS == 1 ) {
         for (J = 1; J <= N - 1; J++) { // 70
            dspr2(CUPLO, N, -E( J ), U( 1, J ), 1, U( 1, J+1 ), 1, WORK );
         } // 70
      }
      WNORM = DLANSP( '1', CUPLO, N, WORK, WORK( LW+1 ) )

      if ( ANORM.GT.WNORM ) {
         RESULT( 1 ) = ( WNORM / ANORM ) / ( N*ULP )
      } else {
         if ( ANORM.LT.ONE ) {
            RESULT( 1 ) = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP )
         } else {
            RESULT( 1 ) = MIN( WNORM / ANORM, DBLE( N ) ) / ( N*ULP )
         }
      }

      // Do Test 2

      // Compute  U U**T - I

      dgemm('N', 'C', N, N, N, ONE, U, LDU, U, LDU, ZERO, WORK, N );

      for (J = 1; J <= N; J++) { // 80
         WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - ONE
      } // 80

      RESULT( 2 ) = MIN( DLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) ), DBLE( N ) ) / ( N*ULP )

      RETURN

      // End of DSBT21

      }
