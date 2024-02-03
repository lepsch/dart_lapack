      SUBROUTINE ZLA_GEAMV( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             ALPHA, BETA;
      int                INCX, INCY, LDA, M, N;
      int                TRANS;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), X( * )
      double             Y( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE, ZERO
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               SYMB_ZERO;
      double             TEMP, SAFE1;
      int                I, INFO, IY, J, JX, KX, KY, LENX, LENY;
      COMPLEX*16         CDUM
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, DLAMCH
      double             DLAMCH;
      // ..
      // .. External Functions ..
      // EXTERNAL ILATRANS
      int                ILATRANS;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, ABS, REAL, DIMAG, SIGN
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function Definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      if ( .NOT.( ( TRANS == ILATRANS( 'N' ) ) || ( TRANS == ILATRANS( 'T' ) ) || ( TRANS == ILATRANS( 'C' ) ) ) ) {
         INFO = 1
      } else if ( M < 0 ) {
         INFO = 2
      } else if ( N < 0 ) {
         INFO = 3
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = 6
      } else if ( INCX == 0 ) {
         INFO = 8
      } else if ( INCY == 0 ) {
         INFO = 11
      }
      if ( INFO != 0 ) {
         xerbla('ZLA_GEAMV ', INFO );
         RETURN
      }

      // Quick return if possible.

      IF( ( M == 0 ) || ( N == 0 ) || ( ( ALPHA == ZERO ) && ( BETA == ONE ) ) ) RETURN

      // Set  LENX  and  LENY, the lengths of the vectors x and y, and set
      // up the start points in  X  and  Y.

      if ( TRANS == ILATRANS( 'N' ) ) {
         LENX = N
         LENY = M
      } else {
         LENX = M
         LENY = N
      }
      if ( INCX > 0 ) {
         KX = 1
      } else {
         KX = 1 - ( LENX - 1 )*INCX
      }
      if ( INCY > 0 ) {
         KY = 1
      } else {
         KY = 1 - ( LENY - 1 )*INCY
      }

      // Set SAFE1 essentially to be the underflow threshold times the
      // number of additions in each row.

      SAFE1 = DLAMCH( 'Safe minimum' )
      SAFE1 = (N+1)*SAFE1

      // Form  y := alpha*abs(A)*abs(x) + beta*abs(y).

      // The O(M*N) SYMB_ZERO tests could be replaced by O(N) queries to
      // the inexact flag.  Still doesn't help change the iteration order
      // to per-column.

      IY = KY
      if ( INCX == 1 ) {
         if ( TRANS == ILATRANS( 'N' ) ) {
            for (I = 1; I <= LENY; I++) {
               if ( BETA == 0.0D+0 ) {
                  SYMB_ZERO = true;
                  Y( IY ) = 0.0D+0
               } else if ( Y( IY ) == 0.0D+0 ) {
                  SYMB_ZERO = true;
               } else {
                  SYMB_ZERO = false;
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               if ( ALPHA != 0.0D+0 ) {
                  for (J = 1; J <= LENX; J++) {
                     TEMP = CABS1( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP
                  }
               }
                if (.NOT.SYMB_ZERO) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) );

               IY = IY + INCY
            }
         } else {
            for (I = 1; I <= LENY; I++) {
               if ( BETA == 0.0D+0 ) {
                  SYMB_ZERO = true;
                  Y( IY ) = 0.0D+0
               } else if ( Y( IY ) == 0.0D+0 ) {
                  SYMB_ZERO = true;
               } else {
                  SYMB_ZERO = false;
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               if ( ALPHA != 0.0D+0 ) {
                  for (J = 1; J <= LENX; J++) {
                     TEMP = CABS1( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO && ( X( J ) == ZERO || TEMP == ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP
                  }
               }
                if (.NOT.SYMB_ZERO) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) );

               IY = IY + INCY
            }
         }
      } else {
         if ( TRANS == ILATRANS( 'N' ) ) {
            for (I = 1; I <= LENY; I++) {
               if ( BETA == 0.0D+0 ) {
                  SYMB_ZERO = true;
                  Y( IY ) = 0.0D+0
               } else if ( Y( IY ) == 0.0D+0 ) {
                  SYMB_ZERO = true;
               } else {
                  SYMB_ZERO = false;
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               if ( ALPHA != 0.0D+0 ) {
                  JX = KX
                  for (J = 1; J <= LENX; J++) {
                     TEMP = CABS1( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO && ( X( JX ) == ZERO || TEMP == ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP
                     JX = JX + INCX
                  }
               }
                if (.NOT.SYMB_ZERO) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) );

               IY = IY + INCY
            }
         } else {
            for (I = 1; I <= LENY; I++) {
               if ( BETA == 0.0D+0 ) {
                  SYMB_ZERO = true;
                  Y( IY ) = 0.0D+0
               } else if ( Y( IY ) == 0.0D+0 ) {
                  SYMB_ZERO = true;
               } else {
                  SYMB_ZERO = false;
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               if ( ALPHA != 0.0D+0 ) {
                  JX = KX
                  for (J = 1; J <= LENX; J++) {
                     TEMP = CABS1( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO && ( X( JX ) == ZERO || TEMP == ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP
                     JX = JX + INCX
                  }
               }
                if (.NOT.SYMB_ZERO) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) );

               IY = IY + INCY
            }
         }

      }

      RETURN

      // End of ZLA_GEAMV

      }
