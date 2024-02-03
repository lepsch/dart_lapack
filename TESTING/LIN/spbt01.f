      SUBROUTINE SPBT01( UPLO, N, KD, A, LDA, AFAC, LDAFAC, RWORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            KD, LDA, LDAFAC, N
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), AFAC( LDAFAC, * ), RWORK( * )
*     ..
*
*  =====================================================================
*
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, K, KC, KLEN, ML, MU
      REAL               ANORM, EPS, T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SDOT, SLAMCH, SLANSB
      EXTERNAL           LSAME, SDOT, SLAMCH, SLANSB
*     ..
*     .. External Subroutines ..
      EXTERNAL           SSCAL, SSYR, STRMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
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
      ANORM = SLANSB( '1', UPLO, N, KD, A, LDA, RWORK )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
*
*     Compute the product U'*U, overwriting U.
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 10 K = N, 1, -1
            KC = MAX( 1, KD+2-K )
            KLEN = KD + 1 - KC
*
*           Compute the (K,K) element of the result.
*
            T = SDOT( KLEN+1, AFAC( KC, K ), 1, AFAC( KC, K ), 1 )
            AFAC( KD+1, K ) = T
*
*           Compute the rest of column K.
*
            IF( KLEN.GT.0 ) CALL STRMV( 'Upper', 'Transpose', 'Non-unit', KLEN, AFAC( KD+1, K-KLEN ), LDAFAC-1, AFAC( KC, K ), 1 )
*
   10    CONTINUE
*
*     UPLO = 'L':  Compute the product L*L', overwriting L.
*
      ELSE
         DO 20 K = N, 1, -1
            KLEN = MIN( KD, N-K )
*
*           Add a multiple of column K of the factor L to each of
*           columns K+1 through N.
*
            IF( KLEN.GT.0 ) CALL SSYR( 'Lower', KLEN, ONE, AFAC( 2, K ), 1, AFAC( 1, K+1 ), LDAFAC-1 )
*
*           Scale column K by the diagonal element.
*
            T = AFAC( 1, K )
            CALL SSCAL( KLEN+1, T, AFAC( 1, K ), 1 )
*
   20    CONTINUE
      END IF
*
*     Compute the difference  L*L' - A  or  U'*U - A.
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 40 J = 1, N
            MU = MAX( 1, KD+2-J )
            DO 30 I = MU, KD + 1
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE
         DO 60 J = 1, N
            ML = MIN( KD+1, N-J+1 )
            DO 50 I = 1, ML
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
*
*     Compute norm( L*L' - A ) / ( N * norm(A) * EPS )
*
      RESID = SLANSB( 'I', UPLO, N, KD, AFAC, LDAFAC, RWORK )
*
      RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS
*
      RETURN
*
*     End of SPBT01
*
      END
