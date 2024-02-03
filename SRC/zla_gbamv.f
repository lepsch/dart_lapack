      SUBROUTINE ZLA_GBAMV( TRANS, M, N, KL, KU, ALPHA, AB, LDAB, X, INCX, BETA, Y, INCY )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             ALPHA, BETA;
      int                INCX, INCY, LDAB, M, N, KL, KU, TRANS;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         AB( LDAB, * ), X( * )
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
      int                I, INFO, IY, J, JX, KX, KY, LENX, LENY, KD, KE;
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
      // .. Statement Functions
      double             CABS1;
      // ..
      // .. Statement Function Definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      if ( .NOT.( ( TRANS.EQ.ILATRANS( 'N' ) ) .OR. ( TRANS.EQ.ILATRANS( 'T' ) ) .OR. ( TRANS.EQ.ILATRANS( 'C' ) ) ) ) {
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
      } else if ( INCX.EQ.0 ) {
         INFO = 8
      } else if ( INCY.EQ.0 ) {
         INFO = 11
      }
      if ( INFO.NE.0 ) {
         xerbla('ZLA_GBAMV ', INFO );
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

      SAFE1 = DLAMCH( 'Safe minimum' )
      SAFE1 = (N+1)*SAFE1

      // Form  y := alpha*abs(A)*abs(x) + beta*abs(y).

      // The O(M*N) SYMB_ZERO tests could be replaced by O(N) queries to
      // the inexact flag.  Still doesn't help change the iteration order
      // to per-column.

      KD = KU + 1
      KE = KL + 1
      IY = KY
      if ( INCX.EQ.1 ) {
         if ( TRANS.EQ.ILATRANS( 'N' ) ) {
            for (I = 1; I <= LENY; I++) {
               if ( BETA .EQ. 0.0D+0 ) {
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0D+0
               } else if ( Y( IY ) .EQ. 0.0D+0 ) {
                  SYMB_ZERO = .TRUE.
               } else {
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               if ( ALPHA .NE. 0.0D+0 ) {
                  DO J = MAX( I-KL, 1 ), MIN( I+KU, LENX )
                     TEMP = CABS1( AB( KD+I-J, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP
                  END DO
               }
                IF ( .NOT.SYMB_ZERO) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         } else {
            for (I = 1; I <= LENY; I++) {
               if ( BETA .EQ. 0.0D+0 ) {
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0D+0
               } else if ( Y( IY ) .EQ. 0.0D+0 ) {
                  SYMB_ZERO = .TRUE.
               } else {
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               if ( ALPHA .NE. 0.0D+0 ) {
                  DO J = MAX( I-KL, 1 ), MIN( I+KU, LENX )
                     TEMP = CABS1( AB( KE-I+J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP
                  END DO
               }
                IF ( .NOT.SYMB_ZERO) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         }
      } else {
         if ( TRANS.EQ.ILATRANS( 'N' ) ) {
            for (I = 1; I <= LENY; I++) {
               if ( BETA .EQ. 0.0D+0 ) {
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0D+0
               } else if ( Y( IY ) .EQ. 0.0D+0 ) {
                  SYMB_ZERO = .TRUE.
               } else {
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               if ( ALPHA .NE. 0.0D+0 ) {
                  JX = KX
                  DO J = MAX( I-KL, 1 ), MIN( I+KU, LENX )
                     TEMP = CABS1( AB( KD+I-J, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( JX ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
               }
                IF ( .NOT.SYMB_ZERO ) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         } else {
            for (I = 1; I <= LENY; I++) {
               if ( BETA .EQ. 0.0D+0 ) {
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0D+0
               } else if ( Y( IY ) .EQ. 0.0D+0 ) {
                  SYMB_ZERO = .TRUE.
               } else {
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               if ( ALPHA .NE. 0.0D+0 ) {
                  JX = KX
                  DO J = MAX( I-KL, 1 ), MIN( I+KU, LENX )
                     TEMP = CABS1( AB( KE-I+J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( JX ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
               }
                IF ( .NOT.SYMB_ZERO ) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         }

      }

      RETURN

      // End of ZLA_GBAMV

      }
