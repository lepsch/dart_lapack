      SUBROUTINE ZGET01( M, N, A, LDA, AFAC, LDAFAC, IPIV, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDAFAC, M, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), AFAC( LDAFAC, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         CONE
      const              CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, J, K;
      double             ANORM, EPS;
      COMPLEX*16         T
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANGE;
      COMPLEX*16         ZDOTU
      // EXTERNAL DLAMCH, ZLANGE, ZDOTU
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMV, ZLASWP, ZSCAL, ZTRMV
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
      ANORM = ZLANGE( '1', M, N, A, LDA, RWORK )

      // Compute the product L*U and overwrite AFAC with the result.
      // A column at a time of the product is obtained, starting with
      // column N.

      DO 10 K = N, 1, -1
         if ( K.GT.M ) {
            ztrmv('Lower', 'No transpose', 'Unit', M, AFAC, LDAFAC, AFAC( 1, K ), 1 );
         } else {

            // Compute elements (K+1:M,K)

            T = AFAC( K, K )
            if ( K+1.LE.M ) {
               zscal(M-K, T, AFAC( K+1, K ), 1 );
               zgemv('No transpose', M-K, K-1, CONE, AFAC( K+1, 1 ), LDAFAC, AFAC( 1, K ), 1, CONE, AFAC( K+1, K ), 1 );
            }

            // Compute the (K,K) element

            AFAC( K, K ) = T + ZDOTU( K-1, AFAC( K, 1 ), LDAFAC, AFAC( 1, K ), 1 )

            // Compute elements (1:K-1,K)

            ztrmv('Lower', 'No transpose', 'Unit', K-1, AFAC, LDAFAC, AFAC( 1, K ), 1 );
         }
   10 CONTINUE
      zlaswp(N, AFAC, LDAFAC, 1, MIN( M, N ), IPIV, -1 );

      // Compute the difference  L*U - A  and store in AFAC.

      for (J = 1; J <= N; J++) { // 30
         for (I = 1; I <= M; I++) { // 20
            AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   20    CONTINUE
   30 CONTINUE

      // Compute norm( L*U - A ) / ( N * norm(A) * EPS )

      RESID = ZLANGE( '1', M, N, AFAC, LDAFAC, RWORK )

      if ( ANORM.LE.ZERO ) {
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      } else {
         RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS
      }

      RETURN

      // End of ZGET01

      }
