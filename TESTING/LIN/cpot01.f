      SUBROUTINE CPOT01( UPLO, N, A, LDA, AFAC, LDAFAC, RWORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDAFAC, N
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), AFAC( LDAFAC, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, K
      REAL               ANORM, EPS, TR
      COMPLEX            TC
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               CLANHE, SLAMCH
      COMPLEX            CDOTC
      EXTERNAL           LSAME, CLANHE, SLAMCH, CDOTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHER, CSCAL, CTRMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          AIMAG, REAL
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0.
*
      IF( N.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
*     Exit with RESID = 1/EPS if ANORM = 0.
*
      EPS = SLAMCH( 'Epsilon' )
      ANORM = CLANHE( '1', UPLO, N, A, LDA, RWORK )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
*
*     Check the imaginary parts of the diagonal elements and return with
*     an error code if any are nonzero.
*
      DO 10 J = 1, N
         IF( AIMAG( AFAC( J, J ) ).NE.ZERO ) THEN
            RESID = ONE / EPS
            RETURN
         END IF
   10 CONTINUE
*
*     Compute the product U**H * U, overwriting U.
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 K = N, 1, -1
*
*           Compute the (K,K) element of the result.
*
            TR = REAL( CDOTC( K, AFAC( 1, K ), 1, AFAC( 1, K ), 1 ) )
            AFAC( K, K ) = TR
*
*           Compute the rest of column K.
*
            CALL CTRMV( 'Upper', 'Conjugate', 'Non-unit', K-1, AFAC, LDAFAC, AFAC( 1, K ), 1 )
*
   20    CONTINUE
*
*     Compute the product L * L**H, overwriting L.
*
      ELSE
         DO 30 K = N, 1, -1
*
*           Add a multiple of column K of the factor L to each of
*           columns K+1 through N.
*
            IF( K+1.LE.N ) CALL CHER( 'Lower', N-K, ONE, AFAC( K+1, K ), 1, AFAC( K+1, K+1 ), LDAFAC )
*
*           Scale column K by the diagonal element.
*
            TC = AFAC( K, K )
            CALL CSCAL( N-K+1, TC, AFAC( K, K ), 1 )
*
   30    CONTINUE
      END IF
*
*     Compute the difference L * L**H - A (or U**H * U - A).
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 50 J = 1, N
            DO 40 I = 1, J - 1
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   40       CONTINUE
            AFAC( J, J ) = AFAC( J, J ) - REAL( A( J, J ) )
   50    CONTINUE
      ELSE
         DO 70 J = 1, N
            AFAC( J, J ) = AFAC( J, J ) - REAL( A( J, J ) )
            DO 60 I = J + 1, N
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   60       CONTINUE
   70    CONTINUE
      END IF
*
*     Compute norm(L*U - A) / ( N * norm(A) * EPS )
*
      RESID = CLANHE( '1', UPLO, N, AFAC, LDAFAC, RWORK )
*
      RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS
*
      RETURN
*
*     End of CPOT01
*
      END
