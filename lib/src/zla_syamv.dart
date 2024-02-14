      void zla_syamv(final int UPLO, final int N, final int ALPHA, final Matrix<double> A_, final int LDA, final int X, final int INCX, final int BETA, final int Y, final int INCY,) {
  final A = A_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double             ALPHA, BETA;
      int                INCX, INCY, LDA, N;
      int                UPLO;
      Complex         A( LDA, * ), X( * );
      double             Y( * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               SYMB_ZERO;
      double             TEMP, SAFE1;
      int                I, INFO, IY, J, JX, KX, KY;
      Complex         ZDUM;
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
      // INTRINSIC MAX, ABS, SIGN, REAL, DIMAG
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function Definitions ..
      CABS1[ZDUM] = ( ZDUM.toDouble() ).abs() + ( DIMAG ( ZDUM ) ).abs();

      // Test the input parameters.

      INFO = 0;
      if ( UPLO != ILAUPLO( 'U' ) && UPLO != ILAUPLO( 'L' ) ) {
         INFO = 1;
      } else if ( N < 0 ) {
         INFO = 2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = 5;
      } else if ( INCX == 0 ) {
         INFO = 7;
      } else if ( INCY == 0 ) {
         INFO = 10;
      }
      if ( INFO != 0 ) {
         xerbla('ZLA_SYAMV', INFO );
         return;
      }

      // Quick return if possible.

      if( ( N == 0 ) || ( ( ALPHA == ZERO ) && ( BETA == ONE ) ) ) return;

      // Set up the start points in  X  and  Y.

      if ( INCX > 0 ) {
         KX = 1;
      } else {
         KX = 1 - ( N - 1 )*INCX;
      }
      if ( INCY > 0 ) {
         KY = 1;
      } else {
         KY = 1 - ( N - 1 )*INCY;
      }

      // Set SAFE1 essentially to be the underflow threshold times the
      // number of additions in each row.

      SAFE1 = dlamch( 'Safe minimum' );
      SAFE1 = (N+1)*SAFE1;

      // Form  y := alpha*abs(A)*abs(x) + beta*abs(y).

      // The O(N^2) SYMB_ZERO tests could be replaced by O(N) queries to
      // the inexact flag.  Still doesn't help change the iteration order
      // to per-column.

      IY = KY;
      if ( INCX == 1 ) {
         if ( UPLO == ILAUPLO( 'U' ) ) {
            for (I = 1; I <= N; I++) {
               if ( BETA == ZERO ) {
                  SYMB_ZERO = true;
                  Y[IY] = 0.0;
               } else if ( Y( IY ) == ZERO ) {
                  SYMB_ZERO = true;
               } else {
                  SYMB_ZERO = false;
                  Y[IY] = BETA * ( Y( IY ) ).abs();
               }
               if ( ALPHA != ZERO ) {
                  for (J = 1; J <= I; J++) {
                     TEMP = CABS1( A( J, I ) );
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO );

                     Y[IY] = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP;
                  }
                  for (J = I+1; J <= N; J++) {
                     TEMP = CABS1( A( I, J ) );
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO );

                     Y[IY] = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP;
                  }
               }
                if ( !SYMB_ZERO) Y( IY ) = Y( IY ) + sign( SAFE1, Y( IY ) );

               IY = IY + INCY;
            }
         } else {
            for (I = 1; I <= N; I++) {
               if ( BETA == ZERO ) {
                  SYMB_ZERO = true;
                  Y[IY] = 0.0;
               } else if ( Y( IY ) == ZERO ) {
                  SYMB_ZERO = true;
               } else {
                  SYMB_ZERO = false;
                  Y[IY] = BETA * ( Y( IY ) ).abs();
               }
               if ( ALPHA != ZERO ) {
                  for (J = 1; J <= I; J++) {
                     TEMP = CABS1( A( I, J ) );
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO );

                     Y[IY] = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP;
                  }
                  for (J = I+1; J <= N; J++) {
                     TEMP = CABS1( A( J, I ) );
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO );

                     Y[IY] = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP;
                  }
               }
                if ( !SYMB_ZERO) Y( IY ) = Y( IY ) + sign( SAFE1, Y( IY ) );

               IY = IY + INCY;
            }
         }
      } else {
         if ( UPLO == ILAUPLO( 'U' ) ) {
            for (I = 1; I <= N; I++) {
               if ( BETA == ZERO ) {
                  SYMB_ZERO = true;
                  Y[IY] = 0.0;
               } else if ( Y( IY ) == ZERO ) {
                  SYMB_ZERO = true;
               } else {
                  SYMB_ZERO = false;
                  Y[IY] = BETA * ( Y( IY ) ).abs();
               }
               JX = KX;
               if ( ALPHA != ZERO ) {
                  for (J = 1; J <= I; J++) {
                     TEMP = CABS1( A( J, I ) );
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO );

                     Y[IY] = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP;
                     JX = JX + INCX;
                  }
                  for (J = I+1; J <= N; J++) {
                     TEMP = CABS1( A( I, J ) );
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO );

                     Y[IY] = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP;
                     JX = JX + INCX;
                  }
               }
                if ( !SYMB_ZERO) Y( IY ) = Y( IY ) + sign( SAFE1, Y( IY ) );

               IY = IY + INCY;
            }
         } else {
            for (I = 1; I <= N; I++) {
               if ( BETA == ZERO ) {
                  SYMB_ZERO = true;
                  Y[IY] = 0.0;
               } else if ( Y( IY ) == ZERO ) {
                  SYMB_ZERO = true;
               } else {
                  SYMB_ZERO = false;
                  Y[IY] = BETA * ( Y( IY ) ).abs();
               }
               JX = KX;
               if ( ALPHA != ZERO ) {
                  for (J = 1; J <= I; J++) {
                     TEMP = CABS1( A( I, J ) );
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO );

                     Y[IY] = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP;
                     JX = JX + INCX;
                  }
                  for (J = I+1; J <= N; J++) {
                     TEMP = CABS1( A( J, I ) );
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO );

                     Y[IY] = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP;
                     JX = JX + INCX;
                  }
               }
                if ( !SYMB_ZERO) Y( IY ) = Y( IY ) + sign( SAFE1, Y( IY ) );

               IY = IY + INCY;
            }
         }

      }

      }
