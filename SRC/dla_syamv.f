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
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
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
      if ( UPLO.NE.ILAUPLO( 'U' ) .AND. UPLO.NE.ILAUPLO( 'L' ) ) {
         INFO = 1
      } else if ( N.LT.0 ) {
         INFO = 2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = 5
      } else if ( INCX.EQ.0 ) {
         INFO = 7
      } else if ( INCY.EQ.0 ) {
         INFO = 10
      }
      if ( INFO.NE.0 ) {
         xerbla('DLA_SYAMV', INFO );
         RETURN
      }

      // Quick return if possible.

      IF( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) ) RETURN

      // Set up the start points in  X  and  Y.

      if ( INCX.GT.0 ) {
         KX = 1
      } else {
         KX = 1 - ( N - 1 )*INCX
      }
      if ( INCY.GT.0 ) {
         KY = 1
      } else {
         KY = 1 - ( N - 1 )*INCY
      }

      // Set SAFE1 essentially to be the underflow threshold times the
      // number of additions in each row.

      SAFE1 = DLAMCH( 'Safe minimum' )
      SAFE1 = (N+1)*SAFE1

      // Form  y := alpha*abs(A)*abs(x) + beta*abs(y).

      // The O(N^2) SYMB_ZERO tests could be replaced by O(N) queries to
      // the inexact flag.  Still doesn't help change the iteration order
      // to per-column.

      IY = KY
      if ( INCX.EQ.1 ) {
         if ( UPLO .EQ. ILAUPLO( 'U' ) ) {
            DO I = 1, N
               if ( BETA .EQ. ZERO ) {
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0D+0
               } else if ( Y( IY ) .EQ. ZERO ) {
                  SYMB_ZERO = .TRUE.
               } else {
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               if ( ALPHA .NE. ZERO ) {
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
               }
                IF ( .NOT.SYMB_ZERO ) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         } else {
            DO I = 1, N
               if ( BETA .EQ. ZERO ) {
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0D+0
               } else if ( Y( IY ) .EQ. ZERO ) {
                  SYMB_ZERO = .TRUE.
               } else {
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               if ( ALPHA .NE. ZERO ) {
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
               }
                IF ( .NOT.SYMB_ZERO ) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         }
      } else {
         if ( UPLO .EQ. ILAUPLO( 'U' ) ) {
            DO I = 1, N
               if ( BETA .EQ. ZERO ) {
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0D+0
               } else if ( Y( IY ) .EQ. ZERO ) {
                  SYMB_ZERO = .TRUE.
               } else {
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               JX = KX
               if ( ALPHA .NE. ZERO ) {
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
               }
                IF ( .NOT.SYMB_ZERO ) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         } else {
            DO I = 1, N
               if ( BETA .EQ. ZERO ) {
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0D+0
               } else if ( Y( IY ) .EQ. ZERO ) {
                  SYMB_ZERO = .TRUE.
               } else {
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               JX = KX
               if ( ALPHA .NE. ZERO ) {
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
               }
                IF ( .NOT.SYMB_ZERO ) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         }

      }

      RETURN

      // End of DLA_SYAMV

      }
