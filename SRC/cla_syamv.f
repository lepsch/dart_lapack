      SUBROUTINE CLA_SYAMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               ALPHA, BETA
      int                INCX, INCY, LDA, N;
      int                UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), X( * )
      REAL               Y( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               SYMB_ZERO;
      REAL               TEMP, SAFE1
      int                I, INFO, IY, J, JX, KX, KY;
      COMPLEX            ZDUM
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, SLAMCH
      REAL               SLAMCH
      // ..
      // .. External Functions ..
      // EXTERNAL ILAUPLO
      int                ILAUPLO;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, ABS, SIGN, REAL, AIMAG
      // ..
      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( REAL ( ZDUM ) ) + ABS( AIMAG ( ZDUM ) )
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
         xerbla('CLA_SYAMV', INFO );
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

      SAFE1 = SLAMCH( 'Safe minimum' )
      SAFE1 = (N+1)*SAFE1

      // Form  y := alpha*abs(A)*abs(x) + beta*abs(y).

      // The O(N^2) SYMB_ZERO tests could be replaced by O(N) queries to
      // the inexact flag.  Still doesn't help change the iteration order
      // to per-column.

      IY = KY
      if ( INCX.EQ.1 ) {
         if ( UPLO .EQ. ILAUPLO( 'U' ) ) {
            for (I = 1; I <= N; I++) {
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
                  for (J = 1; J <= I; J++) {
                     TEMP = CABS1( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP
                  END DO
                  for (J = I+1; J <= N; J++) {
                     TEMP = CABS1( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP
                  END DO
               }
                IF ( .NOT.SYMB_ZERO ) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         } else {
            for (I = 1; I <= N; I++) {
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
                  for (J = 1; J <= I; J++) {
                     TEMP = CABS1( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP
                  END DO
                  for (J = I+1; J <= N; J++) {
                     TEMP = CABS1( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP
                  END DO
               }
                IF ( .NOT.SYMB_ZERO ) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         }
      } else {
         if ( UPLO .EQ. ILAUPLO( 'U' ) ) {
            for (I = 1; I <= N; I++) {
               if ( BETA .EQ. ZERO ) {
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0
               } else if ( Y( IY ) .EQ. ZERO ) {
                  SYMB_ZERO = .TRUE.
               } else {
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               JX = KX
               if ( ALPHA .NE. ZERO ) {
                  for (J = 1; J <= I; J++) {
                     TEMP = CABS1( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
                  for (J = I+1; J <= N; J++) {
                     TEMP = CABS1( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
               }
                IF ( .NOT.SYMB_ZERO ) Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         } else {
            for (I = 1; I <= N; I++) {
               if ( BETA .EQ. ZERO ) {
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0
               } else if ( Y( IY ) .EQ. ZERO ) {
                  SYMB_ZERO = .TRUE.
               } else {
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               }
               JX = KX
               if ( ALPHA .NE. ZERO ) {
                  for (J = 1; J <= I; J++) {
                     TEMP = CABS1( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
                  for (J = I+1; J <= N; J++) {
                     TEMP = CABS1( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND. ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

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

      // End of CLA_SYAMV

      }
