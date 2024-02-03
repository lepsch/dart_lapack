      SUBROUTINE ZSYEQUB( UPLO, N, A, LDA, S, SCOND, AMAX, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, N;
      double             AMAX, SCOND;
      String             UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), WORK( * )
      double             S( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D0, ZERO = 0.0D0 ;
      int                MAX_ITER;
      const              MAX_ITER = 100 ;
      // ..
      // .. Local Scalars ..
      int                I, J, ITER;
      double             AVG, STD, TOL, C0, C1, C2, T, U, SI, D, BASE, SMIN, SMAX, SMLNUM, BIGNUM, SCALE, SUMSQ;
      bool               UP;
      COMPLEX*16         ZDUM
      // ..
      // .. External Functions ..
      double             DLAMCH;
      bool               LSAME;
      // EXTERNAL DLAMCH, LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLASSQ, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, INT, LOG, MAX, MIN, SQRT
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      IF ( .NOT. ( LSAME( UPLO, 'U' ) .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -1
      ELSE IF ( N .LT. 0 ) THEN
         INFO = -2
      ELSE IF ( LDA .LT. MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF ( INFO .NE. 0 ) THEN
         CALL XERBLA( 'ZSYEQUB', -INFO )
         RETURN
      END IF

      UP = LSAME( UPLO, 'U' )
      AMAX = ZERO

      // Quick return if possible.

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
               S( I ) = MAX( S( I ), CABS1( A( I, J ) ) )
               S( J ) = MAX( S( J ), CABS1( A( I, J ) ) )
               AMAX = MAX( AMAX, CABS1( A( I, J ) ) )
            END DO
            S( J ) = MAX( S( J ), CABS1( A( J, J ) ) )
            AMAX = MAX( AMAX, CABS1( A( J, J ) ) )
         END DO
      ELSE
         DO J = 1, N
            S( J ) = MAX( S( J ), CABS1( A( J, J ) ) )
            AMAX = MAX( AMAX, CABS1( A( J, J ) ) )
            DO I = J+1, N
               S( I ) = MAX( S( I ), CABS1( A( I, J ) ) )
               S( J ) = MAX( S( J ), CABS1( A( I, J ) ) )
               AMAX = MAX( AMAX, CABS1( A( I, J ) ) )
            END DO
         END DO
      END IF
      DO J = 1, N
         S( J ) = 1.0D0 / S( J )
      END DO

      TOL = ONE / SQRT( 2.0D0 * N )

      DO ITER = 1, MAX_ITER
         SCALE = 0.0D0
         SUMSQ = 0.0D0
         // beta = |A|s
         DO I = 1, N
            WORK( I ) = ZERO
         END DO
         IF ( UP ) THEN
            DO J = 1, N
               DO I = 1, J-1
                  WORK( I ) = WORK( I ) + CABS1( A( I, J ) ) * S( J )
                  WORK( J ) = WORK( J ) + CABS1( A( I, J ) ) * S( I )
               END DO
               WORK( J ) = WORK( J ) + CABS1( A( J, J ) ) * S( J )
            END DO
         ELSE
            DO J = 1, N
               WORK( J ) = WORK( J ) + CABS1( A( J, J ) ) * S( J )
               DO I = J+1, N
                  WORK( I ) = WORK( I ) + CABS1( A( I, J ) ) * S( J )
                  WORK( J ) = WORK( J ) + CABS1( A( I, J ) ) * S( I )
               END DO
            END DO
         END IF

         // avg = s^T beta / n
         AVG = 0.0D0
         DO I = 1, N
            AVG = AVG + S( I ) * DBLE( WORK( I ) )
         END DO
         AVG = AVG / N

         STD = 0.0D0
         DO I = N+1, 2*N
            WORK( I ) = S( I-N ) * WORK( I-N ) - AVG
         END DO
         CALL ZLASSQ( N, WORK( N+1 ), 1, SCALE, SUMSQ )
         STD = SCALE * SQRT( SUMSQ / N )

         IF ( STD .LT. TOL * AVG ) GOTO 999

         DO I = 1, N
            T = CABS1( A( I, I ) )
            SI = S( I )
            C2 = ( N-1 ) * T
            C1 = ( N-2 ) * ( DBLE( WORK( I ) ) - T*SI )
            C0 = -(T*SI)*SI + 2 * DBLE( WORK( I ) ) * SI - N*AVG
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
                  T = CABS1( A( J, I ) )
                  U = U + S( J )*T
                  WORK( J ) = WORK( J ) + D*T
               END DO
               DO J = I+1,N
                  T = CABS1( A( I, J ) )
                  U = U + S( J )*T
                  WORK( J ) = WORK( J ) + D*T
               END DO
            ELSE
               DO J = 1, I
                  T = CABS1( A( I, J ) )
                  U = U + S( J )*T
                  WORK( J ) = WORK( J ) + D*T
               END DO
               DO J = I+1,N
                  T = CABS1( A( J, I ) )
                  U = U + S( J )*T
                  WORK( J ) = WORK( J ) + D*T
               END DO
            END IF

            AVG = AVG + ( U + DBLE( WORK( I ) ) ) * D / N
            S( I ) = SI
         END DO
      END DO

 999  CONTINUE

      SMLNUM = DLAMCH( 'SAFEMIN' )
      BIGNUM = ONE / SMLNUM
      SMIN = BIGNUM
      SMAX = ZERO
      T = ONE / SQRT( AVG )
      BASE = DLAMCH( 'B' )
      U = ONE / LOG( BASE )
      DO I = 1, N
         S( I ) = BASE ** INT( U * LOG( S( I ) * T ) )
         SMIN = MIN( SMIN, S( I ) )
         SMAX = MAX( SMAX, S( I ) )
      END DO
      SCOND = MAX( SMIN, SMLNUM ) / MIN( SMAX, BIGNUM )

      }
