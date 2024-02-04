      void sla_syamv(UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double               ALPHA, BETA;
      int                INCX, INCY, LDA, N, UPLO;
      // ..
      // .. Array Arguments ..
      double               A( LDA, * ), X( * ), Y( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               SYMB_ZERO;
      double               TEMP, SAFE1;
      int                I, INFO, IY, J, JX, KX, KY;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, SLAMCH
      double               SLAMCH;
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
         xerbla('SLA_SYAMV', INFO );
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

      SAFE1 = SLAMCH( 'Safe minimum' );
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
                     TEMP = ( A( J, I ) ).abs();
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO );

                     Y[IY] = Y( IY ) + ALPHA*( X( J ) ).abs()*TEMP;
                  }
                  for (J = I+1; J <= N; J++) {
                     TEMP = ( A( I, J ) ).abs();
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO );

                     Y[IY] = Y( IY ) + ALPHA*( X( J ) ).abs()*TEMP;
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
                     TEMP = ( A( I, J ) ).abs();
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO );

                     Y[IY] = Y( IY ) + ALPHA*( X( J ) ).abs()*TEMP;
                  }
                  for (J = I+1; J <= N; J++) {
                     TEMP = ( A( J, I ) ).abs();
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO );

                     Y[IY] = Y( IY ) + ALPHA*( X( J ) ).abs()*TEMP;
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
                     TEMP = ( A( J, I ) ).abs();
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO );

                     Y[IY] = Y( IY ) + ALPHA*( X( JX ) ).abs()*TEMP;
                     JX = JX + INCX;
                  }
                  for (J = I+1; J <= N; J++) {
                     TEMP = ( A( I, J ) ).abs();
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO );

                     Y[IY] = Y( IY ) + ALPHA*( X( JX ) ).abs()*TEMP;
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
                     TEMP = ( A( I, J ) ).abs();
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO );

                     Y[IY] = Y( IY ) + ALPHA*( X( JX ) ).abs()*TEMP;
                     JX = JX + INCX;
                  }
                  for (J = I+1; J <= N; J++) {
                     TEMP = ( A( J, I ) ).abs();
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO );

                     Y[IY] = Y( IY ) + ALPHA*( X( JX ) ).abs()*TEMP;
                     JX = JX + INCX;
                  }
               }
                if ( !SYMB_ZERO) Y( IY ) = Y( IY ) + sign( SAFE1, Y( IY ) );

               IY = IY + INCY;
            }
         }

      }

      return;
      }
