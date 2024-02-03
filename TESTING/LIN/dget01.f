      SUBROUTINE DGET01( M, N, A, LDA, AFAC, LDAFAC, IPIV, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDAFAC, M, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             A( LDA, * ), AFAC( LDAFAC, * ), RWORK( * );
      // ..

*  =====================================================================


      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, K;
      double             ANORM, EPS, T;
      // ..
      // .. External Functions ..
      double             DDOT, DLAMCH, DLANGE;
      // EXTERNAL DDOT, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMV, DLASWP, DSCAL, DTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MIN
      // ..
      // .. Executable Statements ..

      // Quick exit if M = 0 or N = 0.

      if ( M.LE.0 .OR. N.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Determine EPS and the norm of A.

      EPS = DLAMCH( 'Epsilon' )
      ANORM = DLANGE( '1', M, N, A, LDA, RWORK )

      // Compute the product L*U and overwrite AFAC with the result.
      // A column at a time of the product is obtained, starting with
      // column N.

      DO 10 K = N, 1, -1
         if ( K.GT.M ) {
            CALL DTRMV( 'Lower', 'No transpose', 'Unit', M, AFAC, LDAFAC, AFAC( 1, K ), 1 )
         } else {

            // Compute elements (K+1:M,K)

            T = AFAC( K, K )
            if ( K+1.LE.M ) {
               CALL DSCAL( M-K, T, AFAC( K+1, K ), 1 )
               CALL DGEMV( 'No transpose', M-K, K-1, ONE, AFAC( K+1, 1 ), LDAFAC, AFAC( 1, K ), 1, ONE, AFAC( K+1, K ), 1 )
            }

            // Compute the (K,K) element

            AFAC( K, K ) = T + DDOT( K-1, AFAC( K, 1 ), LDAFAC, AFAC( 1, K ), 1 )

            // Compute elements (1:K-1,K)

            CALL DTRMV( 'Lower', 'No transpose', 'Unit', K-1, AFAC, LDAFAC, AFAC( 1, K ), 1 )
         }
   10 CONTINUE
      CALL DLASWP( N, AFAC, LDAFAC, 1, MIN( M, N ), IPIV, -1 )

      // Compute the difference  L*U - A  and store in AFAC.

      DO 30 J = 1, N
         DO 20 I = 1, M
            AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   20    CONTINUE
   30 CONTINUE

      // Compute norm( L*U - A ) / ( N * norm(A) * EPS )

      RESID = DLANGE( '1', M, N, AFAC, LDAFAC, RWORK )

      if ( ANORM.LE.ZERO ) {
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      } else {
         RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS
      }

      RETURN

      // End of DGET01

      }
