      SUBROUTINE SLA_GBAMV( TRANS, M, N, KL, KU, ALPHA, AB, LDAB, X, INCX, BETA, Y, INCY )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               ALPHA, BETA
      int                INCX, INCY, LDAB, M, N, KL, KU, TRANS;
      // ..
      // .. Array Arguments ..
      REAL               AB( LDAB, * ), X( * ), Y( * )
      // ..

*  =====================================================================
      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               SYMB_ZERO;
      REAL               TEMP, SAFE1
      int                I, INFO, IY, J, JX, KX, KY, LENX, LENY, KD, KE;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, SLAMCH
      REAL               SLAMCH
      // ..
      // .. External Functions ..
      // EXTERNAL ILATRANS
      int                ILATRANS;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, ABS, SIGN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      if ( .NOT.( ( TRANS == ILATRANS( 'N' ) ) .OR. ( TRANS == ILATRANS( 'T' ) ) .OR. ( TRANS == ILATRANS( 'C' ) ) ) ) {
         INFO = 1
      } else if ( M.LT.0 ) {
         INFO = 2
      } else if ( N.LT.0 ) {
         INFO = 3
      } else if ( KL.LT.0 .OR. KL.GT.M-1 ) {
         INFO = 4
      } else if ( KU.LT.0 .OR. KU.GT.N-1 ) {
         INFO = 5
      } else if ( LDAB.LT.KL+KU+1 ) {
         INFO = 6
      } else if ( INCX == 0 ) {
         INFO = 8
      } else if ( INCY == 0 ) {
         INFO = 11
      }
      if ( INFO != 0 ) {
         xerbla('SLA_GBAMV ', INFO );
         RETURN
      }

      // Quick return if possible.

      IF( ( M == 0 ).OR.( N == 0 ).OR. ( ( ALPHA == ZERO ).AND.( BETA == ONE ) ) ) RETURN

      // Set  LENX  and  LENY, the lengths of the vectors x and y, and set
      // up the start points in  X  and  Y.

      if ( TRANS == ILATRANS( 'N' ) ) {
         LENX = N
         LENY = M
      } else {
         LENX = M
         LENY = N
      }
      if ( INCX.GT.0 ) {
         KX = 1
      } else {
         KX = 1 - ( LENX - 1 )*INCX
      }
      if ( INCY.GT.0 ) {
         KY = 1
      } else {
         KY = 1 - ( LENY - 1 )*INCY
      }

      // Set SAFE1 essentially to be the underflow threshold times the
      // number of additions in each row.

      SAFE1 = SLAMCH( 'Safe minimum' )
      SAFE1 = (N+1)*SAFE1

      // Form  y := alpha*abs(A)*abs(x) + beta*abs(y).

      // The O(M*N) SYMB_ZERO tests could be replaced by O(N) queries to
      // the inexact flag.  Still doesn't help change the iteration order
      // to per-column.

      KD = KU + 1
      KE = KL + 1
      IY = KY
      if ( INCX == 1 ) {
         if ( TRANS == ILATRANS( 'N' ) ) {
            for (I = 1; I <= LENY; I++) {
               if ( BETA == ZERO ) {
                  SYMB_ZERO = true;
                  Y( IY ) = 0.0
               } else if ( Y( IY ) == ZERO ) {
                  SYMB_ZERO = true;
               } else {
                  SYMB_ZERO = false;
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               if ( ALPHA != ZERO ) {
                  DO J = MAX( I-KL, 1 ), MIN( I+KU, LENX )
                     TEMP = ABS( AB( KD+I-J, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( J ) == ZERO .OR. TEMP == ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( J ) )*TEMP
                  }
               }
                if (.NOT.SYMB_ZERO) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) );
               IY = IY + INCY
            }
         } else {
            for (I = 1; I <= LENY; I++) {
               if ( BETA == ZERO ) {
                  SYMB_ZERO = true;
                  Y( IY ) = 0.0
               } else if ( Y( IY ) == ZERO ) {
                  SYMB_ZERO = true;
               } else {
                  SYMB_ZERO = false;
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               if ( ALPHA != ZERO ) {
                  DO J = MAX( I-KL, 1 ), MIN( I+KU, LENX )
                     TEMP = ABS( AB( KE-I+J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( J ) == ZERO .OR. TEMP == ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( J ) )*TEMP
                  }
               }
                if (.NOT.SYMB_ZERO) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) );
               IY = IY + INCY
            }
         }
      } else {
         if ( TRANS == ILATRANS( 'N' ) ) {
            for (I = 1; I <= LENY; I++) {
               if ( BETA == ZERO ) {
                  SYMB_ZERO = true;
                  Y( IY ) = 0.0
               } else if ( Y( IY ) == ZERO ) {
                  SYMB_ZERO = true;
               } else {
                  SYMB_ZERO = false;
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               if ( ALPHA != ZERO ) {
                  JX = KX
                  DO J = MAX( I-KL, 1 ), MIN( I+KU, LENX )
                     TEMP = ABS( AB( KD+I-J, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( JX ) == ZERO .OR. TEMP == ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( JX ) )*TEMP
                     JX = JX + INCX
                  }
               }
                if (.NOT.SYMB_ZERO) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) );

               IY = IY + INCY
            }
         } else {
            for (I = 1; I <= LENY; I++) {
               if ( BETA == ZERO ) {
                  SYMB_ZERO = true;
                  Y( IY ) = 0.0
               } else if ( Y( IY ) == ZERO ) {
                  SYMB_ZERO = true;
               } else {
                  SYMB_ZERO = false;
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               if ( ALPHA != ZERO ) {
                  JX = KX
                  DO J = MAX( I-KL, 1 ), MIN( I+KU, LENX )
                     TEMP = ABS( AB( KE-I+J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( JX ) == ZERO .OR. TEMP == ZERO )

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

      // End of SLA_GBAMV

      }
