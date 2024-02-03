      SUBROUTINE SLA_GEAMV( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               ALPHA, BETA
      int                INCX, INCY, LDA, M, N, TRANS;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), X( * ), Y( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               SYMB_ZERO;
      REAL               TEMP, SAFE1
      int                I, INFO, IY, J, JX, KX, KY, LENX, LENY;
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
      if ( .NOT.( ( TRANS.EQ.ILATRANS( 'N' ) ) .OR. ( TRANS.EQ.ILATRANS( 'T' ) ) .OR. ( TRANS.EQ.ILATRANS( 'C' )) ) ) {
         INFO = 1
      } else if ( M.LT.0 ) {
         INFO = 2
      } else if ( N.LT.0 ) {
         INFO = 3
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = 6
      } else if ( INCX.EQ.0 ) {
         INFO = 8
      } else if ( INCY.EQ.0 ) {
         INFO = 11
      }
      if ( INFO.NE.0 ) {
         xerbla('SLA_GEAMV ', INFO );
         RETURN
      }

      // Quick return if possible.

      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR. ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) ) RETURN

      // Set  LENX  and  LENY, the lengths of the vectors x and y, and set
      // up the start points in  X  and  Y.

      if ( TRANS.EQ.ILATRANS( 'N' ) ) {
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

      IY = KY
      if ( INCX.EQ.1 ) {
         if ( TRANS.EQ.ILATRANS( 'N' ) ) {
            DO I = 1, LENY
               if ( BETA .EQ. ZERO ) {
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0
               } else if ( Y( IY ) .EQ. ZERO ) {
                  SYMB_ZERO = .TRUE.
               } else {
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               if ( ALPHA .NE. ZERO ) {
                  DO J = 1, LENX
                     TEMP = ABS( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( J ) )*TEMP
                  END DO
               }
                IF ( .NOT.SYMB_ZERO ) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         } else {
            DO I = 1, LENY
               if ( BETA .EQ. ZERO ) {
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0
               } else if ( Y( IY ) .EQ. ZERO ) {
                  SYMB_ZERO = .TRUE.
               } else {
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               if ( ALPHA .NE. ZERO ) {
                  DO J = 1, LENX
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
         if ( TRANS.EQ.ILATRANS( 'N' ) ) {
            DO I = 1, LENY
               if ( BETA .EQ. ZERO ) {
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0
               } else if ( Y( IY ) .EQ. ZERO ) {
                  SYMB_ZERO = .TRUE.
               } else {
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               if ( ALPHA .NE. ZERO ) {
                  JX = KX
                  DO J = 1, LENX
                     TEMP = ABS( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( JX ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
               }
                IF (.NOT.SYMB_ZERO) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         } else {
            DO I = 1, LENY
               if ( BETA .EQ. ZERO ) {
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0
               } else if ( Y( IY ) .EQ. ZERO ) {
                  SYMB_ZERO = .TRUE.
               } else {
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               if ( ALPHA .NE. ZERO ) {
                  JX = KX
                  DO J = 1, LENX
                     TEMP = ABS( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( JX ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
               }
                IF (.NOT.SYMB_ZERO) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         }

      }

      RETURN

      // End of SLA_GEAMV

      }
