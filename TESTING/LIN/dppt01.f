      SUBROUTINE DPPT01( UPLO, N, A, AFAC, RWORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                N;
      double             RESID;
*     ..
*     .. Array Arguments ..
      double             A( * ), AFAC( * ), RWORK( * );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      int                I, K, KC, NPP;
      double             ANORM, EPS, T;
*     ..
*     .. External Functions ..
      bool               LSAME;
      double             DDOT, DLAMCH, DLANSP;
      EXTERNAL           LSAME, DDOT, DLAMCH, DLANSP
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSPR, DTPMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
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
      EPS = DLAMCH( 'Epsilon' )
      ANORM = DLANSP( '1', UPLO, N, A, RWORK )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
*
*     Compute the product U'*U, overwriting U.
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         KC = ( N*( N-1 ) ) / 2 + 1
         DO 10 K = N, 1, -1
*
*           Compute the (K,K) element of the result.
*
            T = DDOT( K, AFAC( KC ), 1, AFAC( KC ), 1 )
            AFAC( KC+K-1 ) = T
*
*           Compute the rest of column K.
*
            IF( K.GT.1 ) THEN
               CALL DTPMV( 'Upper', 'Transpose', 'Non-unit', K-1, AFAC, AFAC( KC ), 1 )
               KC = KC - ( K-1 )
            END IF
   10    CONTINUE
*
*     Compute the product L*L', overwriting L.
*
      ELSE
         KC = ( N*( N+1 ) ) / 2
         DO 20 K = N, 1, -1
*
*           Add a multiple of column K of the factor L to each of
*           columns K+1 through N.
*
            IF( K.LT.N ) CALL DSPR( 'Lower', N-K, ONE, AFAC( KC+1 ), 1, AFAC( KC+N-K+1 ) )
*
*           Scale column K by the diagonal element.
*
            T = AFAC( KC )
            CALL DSCAL( N-K+1, T, AFAC( KC ), 1 )
*
            KC = KC - ( N-K+2 )
   20    CONTINUE
      END IF
*
*     Compute the difference  L*L' - A (or U'*U - A).
*
      NPP = N*( N+1 ) / 2
      DO 30 I = 1, NPP
         AFAC( I ) = AFAC( I ) - A( I )
   30 CONTINUE
*
*     Compute norm( L*U - A ) / ( N * norm(A) * EPS )
*
      RESID = DLANSP( '1', UPLO, N, AFAC, RWORK )
*
      RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS
*
      RETURN
*
*     End of DPPT01
*
      END
