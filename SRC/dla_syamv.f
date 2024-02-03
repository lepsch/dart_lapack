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
      if ( UPLO != ILAUPLO( 'U' ) && UPLO != ILAUPLO( 'L' ) ) {
         INFO = 1
      } else if ( N < 0 ) {
         INFO = 2
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = 5
      } else if ( INCX == 0 ) {
         INFO = 7
      } else if ( INCY == 0 ) {
         INFO = 10
      }
      if ( INFO != 0 ) {
         xerbla('DLA_SYAMV', INFO );
         RETURN
      }

      // Quick return if possible.

      IF( ( N == 0 ) || ( ( ALPHA == ZERO ) && ( BETA == ONE ) ) ) RETURN

      // Set up the start points in  X  and  Y.

      if ( INCX > 0 ) {
         KX = 1
      } else {
         KX = 1 - ( N - 1 )*INCX
      }
      if ( INCY > 0 ) {
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
      if ( INCX == 1 ) {
         if ( UPLO == ILAUPLO( 'U' ) ) {
            for (I = 1; I <= N; I++) {
               if ( BETA == ZERO ) {
                  SYMB_ZERO = true;
                  Y( IY ) = 0.0D+0
               } else if ( Y( IY ) == ZERO ) {
                  SYMB_ZERO = true;
               } else {
                  SYMB_ZERO = false;
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               if ( ALPHA != ZERO ) {
                  for (J = 1; J <= I; J++) {
                     TEMP = ABS( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( J ) )*TEMP
                  }
                  for (J = I+1; J <= N; J++) {
                     TEMP = ABS( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( J ) )*TEMP
                  }
               }
                if (.NOT.SYMB_ZERO) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) );

               IY = IY + INCY
            }
         } else {
            for (I = 1; I <= N; I++) {
               if ( BETA == ZERO ) {
                  SYMB_ZERO = true;
                  Y( IY ) = 0.0D+0
               } else if ( Y( IY ) == ZERO ) {
                  SYMB_ZERO = true;
               } else {
                  SYMB_ZERO = false;
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               if ( ALPHA != ZERO ) {
                  for (J = 1; J <= I; J++) {
                     TEMP = ABS( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( J ) )*TEMP
                  }
                  for (J = I+1; J <= N; J++) {
                     TEMP = ABS( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( J ) )*TEMP
                  }
               }
                if (.NOT.SYMB_ZERO) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) );

               IY = IY + INCY
            }
         }
      } else {
         if ( UPLO == ILAUPLO( 'U' ) ) {
            for (I = 1; I <= N; I++) {
               if ( BETA == ZERO ) {
                  SYMB_ZERO = true;
                  Y( IY ) = 0.0D+0
               } else if ( Y( IY ) == ZERO ) {
                  SYMB_ZERO = true;
               } else {
                  SYMB_ZERO = false;
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               JX = KX
               if ( ALPHA != ZERO ) {
                  for (J = 1; J <= I; J++) {
                     TEMP = ABS( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( JX ) )*TEMP
                     JX = JX + INCX
                  }
                  for (J = I+1; J <= N; J++) {
                     TEMP = ABS( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( JX ) )*TEMP
                     JX = JX + INCX
                  }
               }
                if (.NOT.SYMB_ZERO) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) );

               IY = IY + INCY
            }
         } else {
            for (I = 1; I <= N; I++) {
               if ( BETA == ZERO ) {
                  SYMB_ZERO = true;
                  Y( IY ) = 0.0D+0
               } else if ( Y( IY ) == ZERO ) {
                  SYMB_ZERO = true;
               } else {
                  SYMB_ZERO = false;
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               JX = KX
               if ( ALPHA != ZERO ) {
                  for (J = 1; J <= I; J++) {
                     TEMP = ABS( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( JX ) )*TEMP
                     JX = JX + INCX
                  }
                  for (J = I+1; J <= N; J++) {
                     TEMP = ABS( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( JX ) )*TEMP
                     JX = JX + INCX
                  }
               }
                if (.NOT.SYMB_ZERO) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) );

               IY = IY + INCY
            }
         }

      }

      RETURN

      // End of DLA_SYAMV

      }
