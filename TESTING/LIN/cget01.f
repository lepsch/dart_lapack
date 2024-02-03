      SUBROUTINE CGET01( M, N, A, LDA, AFAC, LDAFAC, IPIV, RWORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                LDA, LDAFAC, M, N;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), AFAC( LDAFAC, * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CONE
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ) )
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
*
      // Quick exit if M = 0 or N = 0.
*
      IF( M.LE.0 .OR. N.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
      // Determine EPS and the norm of A.
*
      EPS = SLAMCH( 'Epsilon' )
      ANORM = CLANGE( '1', M, N, A, LDA, RWORK )
*
      // Compute the product L*U and overwrite AFAC with the result.
      // A column at a time of the product is obtained, starting with
      // column N.
*
      DO 10 K = N, 1, -1
         IF( K.GT.M ) THEN
            CALL CTRMV( 'Lower', 'No transpose', 'Unit', M, AFAC, LDAFAC, AFAC( 1, K ), 1 )
         ELSE
*
            // Compute elements (K+1:M,K)
*
            T = AFAC( K, K )
            IF( K+1.LE.M ) THEN
               CALL CSCAL( M-K, T, AFAC( K+1, K ), 1 )
               CALL CGEMV( 'No transpose', M-K, K-1, CONE, AFAC( K+1, 1 ), LDAFAC, AFAC( 1, K ), 1, CONE, AFAC( K+1, K ), 1 )
            END IF
*
            // Compute the (K,K) element
*
            AFAC( K, K ) = T + CDOTU( K-1, AFAC( K, 1 ), LDAFAC, AFAC( 1, K ), 1 )
*
            // Compute elements (1:K-1,K)
*
            CALL CTRMV( 'Lower', 'No transpose', 'Unit', K-1, AFAC, LDAFAC, AFAC( 1, K ), 1 )
         END IF
   10 CONTINUE
      CALL CLASWP( N, AFAC, LDAFAC, 1, MIN( M, N ), IPIV, -1 )
*
      // Compute the difference  L*U - A  and store in AFAC.
*
      DO 30 J = 1, N
         DO 20 I = 1, M
            AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   20    CONTINUE
   30 CONTINUE
*
      // Compute norm( L*U - A ) / ( N * norm(A) * EPS )
*
      RESID = CLANGE( '1', M, N, AFAC, LDAFAC, RWORK )
*
      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      ELSE
         RESID = ( ( RESID/REAL( N ) )/ANORM ) / EPS
      END IF
*
      RETURN
*
      // End of CGET01
*
      END
