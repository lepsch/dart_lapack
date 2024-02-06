      void sla_gbamv(TRANS, M, N, KL, KU, ALPHA, AB, LDAB, X, INCX, BETA, Y, INCY ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double               ALPHA, BETA;
      int                INCX, INCY, LDAB, M, N, KL, KU, TRANS;
      double               AB( LDAB, * ), X( * ), Y( * );
      // ..

// =====================================================================
      // .. Parameters ..
      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               SYMB_ZERO;
      double               TEMP, SAFE1;
      int                I, INFO, IY, J, JX, KX, KY, LENX, LENY, KD, KE;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, SLAMCH
      double               SLAMCH;
      // ..
      // .. External Functions ..
      // EXTERNAL ILATRANS
      int                ILATRANS;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, ABS, SIGN

      // Test the input parameters.

      INFO = 0;
      if ( !( ( TRANS == ILATRANS( 'N' ) ) || ( TRANS == ILATRANS( 'T' ) ) || ( TRANS == ILATRANS( 'C' ) ) ) ) {
         INFO = 1;
      } else if ( M < 0 ) {
         INFO = 2;
      } else if ( N < 0 ) {
         INFO = 3;
      } else if ( KL < 0 || KL > M-1 ) {
         INFO = 4;
      } else if ( KU < 0 || KU > N-1 ) {
         INFO = 5;
      } else if ( LDAB < KL+KU+1 ) {
         INFO = 6;
      } else if ( INCX == 0 ) {
         INFO = 8;
      } else if ( INCY == 0 ) {
         INFO = 11;
      }
      if ( INFO != 0 ) {
         xerbla('SLA_GBAMV ', INFO );
         return;
      }

      // Quick return if possible.

      if( ( M == 0 ) || ( N == 0 ) || ( ( ALPHA == ZERO ) && ( BETA == ONE ) ) ) return;

      // Set  LENX  and  LENY, the lengths of the vectors x and y, and set
      // up the start points in  X  and  Y.

      if ( TRANS == ILATRANS( 'N' ) ) {
         LENX = N;
         LENY = M;
      } else {
         LENX = M;
         LENY = N;
      }
      if ( INCX > 0 ) {
         KX = 1;
      } else {
         KX = 1 - ( LENX - 1 )*INCX;
      }
      if ( INCY > 0 ) {
         KY = 1;
      } else {
         KY = 1 - ( LENY - 1 )*INCY;
      }

      // Set SAFE1 essentially to be the underflow threshold times the
      // number of additions in each row.

      SAFE1 = SLAMCH( 'Safe minimum' );
      SAFE1 = (N+1)*SAFE1;

      // Form  y := alpha*abs(A)*abs(x) + beta*abs(y).

      // The O(M*N) SYMB_ZERO tests could be replaced by O(N) queries to
      // the inexact flag.  Still doesn't help change the iteration order
      // to per-column.

      KD = KU + 1;
      KE = KL + 1;
      IY = KY;
      if ( INCX == 1 ) {
         if ( TRANS == ILATRANS( 'N' ) ) {
            for (I = 1; I <= LENY; I++) {
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
                  for (J = max( I-KL, 1 ); J <= min( I+KU, LENX ); J++) {
                     TEMP = ( AB( KD+I-J, J ) ).abs();
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO );

                     Y[IY] = Y( IY ) + ALPHA*( X( J ) ).abs()*TEMP;
                  }
               }
                if ( !SYMB_ZERO) Y( IY ) = Y( IY ) + sign( SAFE1, Y( IY ) );
               IY = IY + INCY;
            }
         } else {
            for (I = 1; I <= LENY; I++) {
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
                  for (J = max( I-KL, 1 ); J <= min( I+KU, LENX ); J++) {
                     TEMP = ( AB( KE-I+J, I ) ).abs();
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO );

                     Y[IY] = Y( IY ) + ALPHA*( X( J ) ).abs()*TEMP;
                  }
               }
                if ( !SYMB_ZERO) Y( IY ) = Y( IY ) + sign( SAFE1, Y( IY ) );
               IY = IY + INCY;
            }
         }
      } else {
         if ( TRANS == ILATRANS( 'N' ) ) {
            for (I = 1; I <= LENY; I++) {
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
                  JX = KX;
                  for (J = max( I-KL, 1 ); J <= min( I+KU, LENX ); J++) {
                     TEMP = ( AB( KD+I-J, J ) ).abs();
                     SYMB_ZERO = SYMB_ZERO && ( X( JX ) == ZERO || TEMP == ZERO );

                     Y[IY] = Y( IY ) + ALPHA*( X( JX ) ).abs()*TEMP;
                     JX = JX + INCX;
                  }
               }
                if ( !SYMB_ZERO) Y( IY ) = Y( IY ) + sign( SAFE1, Y( IY ) );

               IY = IY + INCY;
            }
         } else {
            for (I = 1; I <= LENY; I++) {
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
                  JX = KX;
                  for (J = max( I-KL, 1 ); J <= min( I+KU, LENX ); J++) {
                     TEMP = ( AB( KE-I+J, I ) ).abs();
                     SYMB_ZERO = SYMB_ZERO && ( X( JX ) == ZERO || TEMP == ZERO );

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
