      SUBROUTINE ZPOT01( UPLO, N, A, LDA, AFAC, LDAFAC, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDAFAC, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), AFAC( LDAFAC, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      // ..
      // .. Local Scalars ..
      int                I, J, K;
      double             ANORM, EPS, TR;
      COMPLEX*16         TC
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANHE;
      COMPLEX*16         ZDOTC
      // EXTERNAL LSAME, DLAMCH, ZLANHE, ZDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZHER, ZSCAL, ZTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DIMAG
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      IF( N.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANHE( '1', UPLO, N, A, LDA, RWORK )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF

      // Check the imaginary parts of the diagonal elements and return with
      // an error code if any are nonzero.

      DO 10 J = 1, N
         IF( DIMAG( AFAC( J, J ) ).NE.ZERO ) THEN
            RESID = ONE / EPS
            RETURN
         END IF
   10 CONTINUE

      // Compute the product U**H * U, overwriting U.

      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 K = N, 1, -1

            // Compute the (K,K) element of the result.

            TR = DBLE( ZDOTC( K, AFAC( 1, K ), 1, AFAC( 1, K ), 1 ) )
            AFAC( K, K ) = TR

            // Compute the rest of column K.

            CALL ZTRMV( 'Upper', 'Conjugate', 'Non-unit', K-1, AFAC, LDAFAC, AFAC( 1, K ), 1 )

   20    CONTINUE

      // Compute the product L * L**H, overwriting L.

      ELSE
         DO 30 K = N, 1, -1

            // Add a multiple of column K of the factor L to each of
            // columns K+1 through N.

            IF( K+1.LE.N ) CALL ZHER( 'Lower', N-K, ONE, AFAC( K+1, K ), 1, AFAC( K+1, K+1 ), LDAFAC )

            // Scale column K by the diagonal element.

            TC = AFAC( K, K )
            CALL ZSCAL( N-K+1, TC, AFAC( K, K ), 1 )

   30    CONTINUE
      END IF

      // Compute the difference L * L**H - A (or U**H * U - A).

      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 50 J = 1, N
            DO 40 I = 1, J - 1
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   40       CONTINUE
            AFAC( J, J ) = AFAC( J, J ) - DBLE( A( J, J ) )
   50    CONTINUE
      ELSE
         DO 70 J = 1, N
            AFAC( J, J ) = AFAC( J, J ) - DBLE( A( J, J ) )
            DO 60 I = J + 1, N
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   60       CONTINUE
   70    CONTINUE
      END IF

      // Compute norm(L*U - A) / ( N * norm(A) * EPS )

      RESID = ZLANHE( '1', UPLO, N, AFAC, LDAFAC, RWORK )

      RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS

      RETURN

      // End of ZPOT01

      END
