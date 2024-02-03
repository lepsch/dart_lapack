      SUBROUTINE SSYEQUB( UPLO, N, A, LDA, S, SCOND, AMAX, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, LDA, N
      REAL               AMAX, SCOND
      String             UPLO;
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), S( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E0, ZERO = 0.0E0 )
      int                MAX_ITER
      PARAMETER          ( MAX_ITER = 100 )
*     ..
*     .. Local Scalars ..
      int                I, J, ITER
      REAL               AVG, STD, TOL, C0, C1, C2, T, U, SI, D, BASE, SMIN, SMAX, SMLNUM, BIGNUM, SCALE, SUMSQ
      LOGICAL            UP
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      LOGICAL            LSAME
      EXTERNAL           LSAME, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLASSQ, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF ( .NOT. ( LSAME( UPLO, 'U' ) .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -1
      ELSE IF ( N .LT. 0 ) THEN
         INFO = -2
      ELSE IF ( LDA .LT. MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF ( INFO .NE. 0 ) THEN
         CALL XERBLA( 'SSYEQUB', -INFO )
         RETURN
      END IF

      UP = LSAME( UPLO, 'U' )
      AMAX = ZERO
*
*     Quick return if possible.
*
      IF ( N .EQ. 0 ) THEN
         SCOND = ONE
         RETURN
      END IF

      DO I = 1, N
         S( I ) = ZERO
      END DO

      AMAX = ZERO
      IF ( UP ) THEN
         DO J = 1, N
            DO I = 1, J-1
               S( I ) = MAX( S( I ), ABS( A( I, J ) ) )
               S( J ) = MAX( S( J ), ABS( A( I, J ) ) )
               AMAX = MAX( AMAX, ABS( A( I, J ) ) )
            END DO
            S( J ) = MAX( S( J ), ABS( A( J, J ) ) )
            AMAX = MAX( AMAX, ABS( A( J, J ) ) )
         END DO
      ELSE
         DO J = 1, N
            S( J ) = MAX( S( J ), ABS( A( J, J ) ) )
            AMAX = MAX( AMAX, ABS( A( J, J ) ) )
            DO I = J+1, N
               S( I ) = MAX( S( I ), ABS( A( I, J ) ) )
               S( J ) = MAX( S( J ), ABS( A( I, J ) ) )
               AMAX = MAX( AMAX, ABS( A( I, J ) ) )
            END DO
         END DO
      END IF
      DO J = 1, N
         S( J ) = 1.0E0 / S( J )
      END DO

      TOL = ONE / SQRT( 2.0E0 * N )

      DO ITER = 1, MAX_ITER
         SCALE = 0.0E0
         SUMSQ = 0.0E0
*        beta = |A|s
         DO I = 1, N
            WORK( I ) = ZERO
         END DO
         IF ( UP ) THEN
            DO J = 1, N
               DO I = 1, J-1
                  WORK( I ) = WORK( I ) + ABS( A( I, J ) ) * S( J )
                  WORK( J ) = WORK( J ) + ABS( A( I, J ) ) * S( I )
               END DO
               WORK( J ) = WORK( J ) + ABS( A( J, J ) ) * S( J )
            END DO
         ELSE
            DO J = 1, N
               WORK( J ) = WORK( J ) + ABS( A( J, J ) ) * S( J )
               DO I = J+1, N
                  WORK( I ) = WORK( I ) + ABS( A( I, J ) ) * S( J )
                  WORK( J ) = WORK( J ) + ABS( A( I, J ) ) * S( I )
               END DO
            END DO
         END IF

*        avg = s^T beta / n
         AVG = 0.0E0
         DO I = 1, N
            AVG = AVG + S( I )*WORK( I )
         END DO
         AVG = AVG / N

         STD = 0.0E0
         DO I = N+1, 2*N
            WORK( I ) = S( I-N ) * WORK( I-N ) - AVG
         END DO
         CALL SLASSQ( N, WORK( N+1 ), 1, SCALE, SUMSQ )
         STD = SCALE * SQRT( SUMSQ / N )

         IF ( STD .LT. TOL * AVG ) GOTO 999

         DO I = 1, N
            T = ABS( A( I, I ) )
            SI = S( I )
            C2 = ( N-1 ) * T
            C1 = ( N-2 ) * ( WORK( I ) - T*SI )
            C0 = -(T*SI)*SI + 2*WORK( I )*SI - N*AVG
            D = C1*C1 - 4*C0*C2

            IF ( D .LE. 0 ) THEN
               INFO = -1
               RETURN
            END IF
            SI = -2*C0 / ( C1 + SQRT( D ) )

            D = SI - S( I )
            U = ZERO
            IF ( UP ) THEN
               DO J = 1, I
                  T = ABS( A( J, I ) )
                  U = U + S( J )*T
                  WORK( J ) = WORK( J ) + D*T
               END DO
               DO J = I+1,N
                  T = ABS( A( I, J ) )
                  U = U + S( J )*T
                  WORK( J ) = WORK( J ) + D*T
               END DO
            ELSE
               DO J = 1, I
                  T = ABS( A( I, J ) )
                  U = U + S( J )*T
                  WORK( J ) = WORK( J ) + D*T
               END DO
               DO J = I+1,N
                  T = ABS( A( J, I ) )
                  U = U + S( J )*T
                  WORK( J ) = WORK( J ) + D*T
               END DO
            END IF

            AVG = AVG + ( U + WORK( I ) ) * D / N
            S( I ) = SI
         END DO
      END DO

 999  CONTINUE

      SMLNUM = SLAMCH( 'SAFEMIN' )
      BIGNUM = ONE / SMLNUM
      SMIN = BIGNUM
      SMAX = ZERO
      T = ONE / SQRT( AVG )
      BASE = SLAMCH( 'B' )
      U = ONE / LOG( BASE )
      DO I = 1, N
         S( I ) = BASE ** INT( U * LOG( S( I ) * T ) )
         SMIN = MIN( SMIN, S( I ) )
         SMAX = MAX( SMAX, S( I ) )
      END DO
      SCOND = MAX( SMIN, SMLNUM ) / MIN( SMAX, BIGNUM )
*
      END
