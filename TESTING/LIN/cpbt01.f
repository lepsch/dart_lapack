      SUBROUTINE CPBT01( UPLO, N, KD, A, LDA, AFAC, LDAFAC, RWORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                KD, LDA, LDAFAC, N;
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), AFAC( LDAFAC, * )
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
      int                I, J, K, KC, KLEN, ML, MU;
      REAL               AKK, ANORM, EPS
*     ..
*     .. External Functions ..
      bool               LSAME;
      REAL               CLANHB, SLAMCH
      COMPLEX            CDOTC
      EXTERNAL           LSAME, CLANHB, SLAMCH, CDOTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHER, CSSCAL, CTRMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          AIMAG, MAX, MIN, REAL
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
      ANORM = CLANHB( '1', UPLO, N, KD, A, LDA, RWORK )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
*
*     Check the imaginary parts of the diagonal elements and return with
*     an error code if any are nonzero.
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 10 J = 1, N
            IF( AIMAG( AFAC( KD+1, J ) ).NE.ZERO ) THEN
               RESID = ONE / EPS
               RETURN
            END IF
   10    CONTINUE
      ELSE
         DO 20 J = 1, N
            IF( AIMAG( AFAC( 1, J ) ).NE.ZERO ) THEN
               RESID = ONE / EPS
               RETURN
            END IF
   20    CONTINUE
      END IF
*
*     Compute the product U'*U, overwriting U.
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 30 K = N, 1, -1
            KC = MAX( 1, KD+2-K )
            KLEN = KD + 1 - KC
*
*           Compute the (K,K) element of the result.
*
            AKK = REAL( CDOTC( KLEN+1, AFAC( KC, K ), 1, AFAC( KC, K ), 1 ) )
            AFAC( KD+1, K ) = AKK
*
*           Compute the rest of column K.
*
            IF( KLEN.GT.0 ) CALL CTRMV( 'Upper', 'Conjugate', 'Non-unit', KLEN, AFAC( KD+1, K-KLEN ), LDAFAC-1, AFAC( KC, K ), 1 )
*
   30    CONTINUE
*
*     UPLO = 'L':  Compute the product L*L', overwriting L.
*
      ELSE
         DO 40 K = N, 1, -1
            KLEN = MIN( KD, N-K )
*
*           Add a multiple of column K of the factor L to each of
*           columns K+1 through N.
*
            IF( KLEN.GT.0 ) CALL CHER( 'Lower', KLEN, ONE, AFAC( 2, K ), 1, AFAC( 1, K+1 ), LDAFAC-1 )
*
*           Scale column K by the diagonal element.
*
            AKK = REAL( AFAC( 1, K ) )
            CALL CSSCAL( KLEN+1, AKK, AFAC( 1, K ), 1 )
*
   40    CONTINUE
      END IF
*
*     Compute the difference  L*L' - A  or  U'*U - A.
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 60 J = 1, N
            MU = MAX( 1, KD+2-J )
            DO 50 I = MU, KD + 1
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   50       CONTINUE
   60    CONTINUE
      ELSE
         DO 80 J = 1, N
            ML = MIN( KD+1, N-J+1 )
            DO 70 I = 1, ML
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   70       CONTINUE
   80    CONTINUE
      END IF
*
*     Compute norm( L*L' - A ) / ( N * norm(A) * EPS )
*
      RESID = CLANHB( '1', UPLO, N, KD, AFAC, LDAFAC, RWORK )
*
      RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS
*
      RETURN
*
*     End of CPBT01
*
      END
