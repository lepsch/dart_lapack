      SUBROUTINE DLA_SYAMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             ALPHA, BETA;
      int                INCX, INCY, LDA, N, UPLO;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), X( * ), Y( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      // ..
      // .. Local Scalars ..
      bool               SYMB_ZERO;
      double             TEMP, SAFE1;
      int                I, INFO, IY, J, JX, KX, KY;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, DLAMCH
      double             DLAMCH;
      // ..
      // .. External Functions ..
      // EXTERNAL ILAUPLO
      int                ILAUPLO;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, ABS, SIGN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      IF     ( UPLO.NE.ILAUPLO( 'U' ) .AND. UPLO.NE.ILAUPLO( 'L' ) ) THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 5
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 7
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DLA_SYAMV', INFO )
         RETURN
      END IF

      // Quick return if possible.

      IF( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) ) RETURN

      // Set up the start points in  X  and  Y.

      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( N - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 )*INCY
      END IF

      // Set SAFE1 essentially to be the underflow threshold times the
      // number of additions in each row.

      SAFE1 = DLAMCH( 'Safe minimum' )
      SAFE1 = (N+1)*SAFE1

      // Form  y := alpha*abs(A)*abs(x) + beta*abs(y).

      // The O(N^2) SYMB_ZERO tests could be replaced by O(N) queries to
     t // he inexact flag.  Still doesn't help change the iteration order
     t // o per-column.

      IY = KY
      IF ( INCX.EQ.1 ) THEN
         IF ( UPLO .EQ. ILAUPLO( 'U' ) ) THEN
            DO I = 1, N
               IF ( BETA .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0D+0
               ELSE IF ( Y( IY ) .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
               ELSE
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               END IF
               IF ( ALPHA .NE. ZERO ) THEN
                  DO J = 1, I
                     TEMP = ABS( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( J ) )*TEMP
                  END DO
                  DO J = I+1, N
                     TEMP = ABS( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( J ) )*TEMP
                  END DO
               END IF
                IF ( .NOT.SYMB_ZERO ) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         ELSE
            DO I = 1, N
               IF ( BETA .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0D+0
               ELSE IF ( Y( IY ) .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
               ELSE
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               END IF
               IF ( ALPHA .NE. ZERO ) THEN
                  DO J = 1, I
                     TEMP = ABS( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( J ) )*TEMP
                  END DO
                  DO J = I+1, N
                     TEMP = ABS( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( J ) )*TEMP
                  END DO
               END IF
                IF ( .NOT.SYMB_ZERO ) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         END IF
      ELSE
         IF ( UPLO .EQ. ILAUPLO( 'U' ) ) THEN
            DO I = 1, N
               IF ( BETA .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0D+0
               ELSE IF ( Y( IY ) .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
               ELSE
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               END IF
               JX = KX
               IF ( ALPHA .NE. ZERO ) THEN
                  DO J = 1, I
                     TEMP = ABS( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
                  DO J = I+1, N
                     TEMP = ABS( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
               END IF
                IF ( .NOT.SYMB_ZERO ) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         ELSE
            DO I = 1, N
               IF ( BETA .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0D+0
               ELSE IF ( Y( IY ) .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
               ELSE
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               END IF
               JX = KX
               IF ( ALPHA .NE. ZERO ) THEN
                  DO J = 1, I
                     TEMP = ABS( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
                  DO J = I+1, N
                     TEMP = ABS( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
               END IF
                IF ( .NOT.SYMB_ZERO ) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         END IF

      END IF

      RETURN

      // End of DLA_SYAMV

      END
