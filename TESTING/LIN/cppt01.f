      SUBROUTINE CPPT01( UPLO, N, A, AFAC, RWORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                N
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( * ), AFAC( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      int                I, K, KC
      REAL               ANORM, EPS, TR
      COMPLEX            TC
*     ..
*     .. External Functions ..
      bool               LSAME;
      REAL               CLANHP, SLAMCH
      COMPLEX            CDOTC
      EXTERNAL           LSAME, CLANHP, SLAMCH, CDOTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHPR, CSCAL, CTPMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          AIMAG, REAL
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0
*
      IF( N.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
*     Exit with RESID = 1/EPS if ANORM = 0.
*
      EPS = SLAMCH( 'Epsilon' )
      ANORM = CLANHP( '1', UPLO, N, A, RWORK )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
*
*     Check the imaginary parts of the diagonal elements and return with
*     an error code if any are nonzero.
*
      KC = 1
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 10 K = 1, N
            IF( AIMAG( AFAC( KC ) ).NE.ZERO ) THEN
               RESID = ONE / EPS
               RETURN
            END IF
            KC = KC + K + 1
   10    CONTINUE
      ELSE
         DO 20 K = 1, N
            IF( AIMAG( AFAC( KC ) ).NE.ZERO ) THEN
               RESID = ONE / EPS
               RETURN
            END IF
            KC = KC + N - K + 1
   20    CONTINUE
      END IF
*
*     Compute the product U'*U, overwriting U.
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         KC = ( N*( N-1 ) ) / 2 + 1
         DO 30 K = N, 1, -1
*
*           Compute the (K,K) element of the result.
*
            TR = REAL( CDOTC( K, AFAC( KC ), 1, AFAC( KC ), 1 ) )
            AFAC( KC+K-1 ) = TR
*
*           Compute the rest of column K.
*
            IF( K.GT.1 ) THEN
               CALL CTPMV( 'Upper', 'Conjugate', 'Non-unit', K-1, AFAC, AFAC( KC ), 1 )
               KC = KC - ( K-1 )
            END IF
   30    CONTINUE
*
*        Compute the difference  L*L' - A
*
         KC = 1
         DO 50 K = 1, N
            DO 40 I = 1, K - 1
               AFAC( KC+I-1 ) = AFAC( KC+I-1 ) - A( KC+I-1 )
   40       CONTINUE
            AFAC( KC+K-1 ) = AFAC( KC+K-1 ) - REAL( A( KC+K-1 ) )
            KC = KC + K
   50    CONTINUE
*
*     Compute the product L*L', overwriting L.
*
      ELSE
         KC = ( N*( N+1 ) ) / 2
         DO 60 K = N, 1, -1
*
*           Add a multiple of column K of the factor L to each of
*           columns K+1 through N.
*
            IF( K.LT.N ) CALL CHPR( 'Lower', N-K, ONE, AFAC( KC+1 ), 1, AFAC( KC+N-K+1 ) )
*
*           Scale column K by the diagonal element.
*
            TC = AFAC( KC )
            CALL CSCAL( N-K+1, TC, AFAC( KC ), 1 )
*
            KC = KC - ( N-K+2 )
   60    CONTINUE
*
*        Compute the difference  U'*U - A
*
         KC = 1
         DO 80 K = 1, N
            AFAC( KC ) = AFAC( KC ) - REAL( A( KC ) )
            DO 70 I = K + 1, N
               AFAC( KC+I-K ) = AFAC( KC+I-K ) - A( KC+I-K )
   70       CONTINUE
            KC = KC + N - K + 1
   80    CONTINUE
      END IF
*
*     Compute norm( L*U - A ) / ( N * norm(A) * EPS )
*
      RESID = CLANHP( '1', UPLO, N, AFAC, RWORK )
*
      RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS
*
      RETURN
*
*     End of CPPT01
*
      END
