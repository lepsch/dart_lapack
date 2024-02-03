      SUBROUTINE CGET01( M, N, A, LDA, AFAC, LDAFAC, IPIV, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDAFAC, M, N;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), AFAC( LDAFAC, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            CONE
      const              CONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, J, K;
      REAL               ANORM, EPS
      COMPLEX            T
      // ..
      // .. External Functions ..
      REAL               CLANGE, SLAMCH
      COMPLEX            CDOTU
      // EXTERNAL CLANGE, SLAMCH, CDOTU
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV, CLASWP, CSCAL, CTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN, REAL
      // ..
      // .. Executable Statements ..

      // Quick exit if M = 0 or N = 0.

      if ( M.LE.0 .OR. N.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Determine EPS and the norm of A.

      EPS = SLAMCH( 'Epsilon' )
      ANORM = CLANGE( '1', M, N, A, LDA, RWORK )

      // Compute the product L*U and overwrite AFAC with the result.
      // A column at a time of the product is obtained, starting with
      // column N.

      DO 10 K = N, 1, -1
         if ( K.GT.M ) {
            ctrmv('Lower', 'No transpose', 'Unit', M, AFAC, LDAFAC, AFAC( 1, K ), 1 );
         } else {

            // Compute elements (K+1:M,K)

            T = AFAC( K, K )
            if ( K+1.LE.M ) {
               cscal(M-K, T, AFAC( K+1, K ), 1 );
               cgemv('No transpose', M-K, K-1, CONE, AFAC( K+1, 1 ), LDAFAC, AFAC( 1, K ), 1, CONE, AFAC( K+1, K ), 1 );
            }

            // Compute the (K,K) element

            AFAC( K, K ) = T + CDOTU( K-1, AFAC( K, 1 ), LDAFAC, AFAC( 1, K ), 1 )

            // Compute elements (1:K-1,K)

            ctrmv('Lower', 'No transpose', 'Unit', K-1, AFAC, LDAFAC, AFAC( 1, K ), 1 );
         }
   10 CONTINUE
      claswp(N, AFAC, LDAFAC, 1, MIN( M, N ), IPIV, -1 );

      // Compute the difference  L*U - A  and store in AFAC.

      for (J = 1; J <= N; J++) { // 30
         for (I = 1; I <= M; I++) { // 20
            AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   20    CONTINUE
   30 CONTINUE

      // Compute norm( L*U - A ) / ( N * norm(A) * EPS )

      RESID = CLANGE( '1', M, N, AFAC, LDAFAC, RWORK )

      if ( ANORM.LE.ZERO ) {
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      } else {
         RESID = ( ( RESID/REAL( N ) )/ANORM ) / EPS
      }

      RETURN

      // End of CGET01

      }
