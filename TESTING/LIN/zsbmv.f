      SUBROUTINE ZSBMV( UPLO, N, K, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INCX, INCY, K, LDA, N;
      COMPLEX*16         ALPHA, BETA
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), X( * ), Y( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE
      const              ONE = ( 1.0D+0, 0.0D+0 ) ;
      COMPLEX*16         ZERO
      const              ZERO = ( 0.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, IX, IY, J, JX, JY, KPLUS1, KX, KY, L;
      COMPLEX*16         TEMP1, TEMP2
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      if ( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = 1
      } else if ( N.LT.0 ) {
         INFO = 2
      } else if ( K.LT.0 ) {
         INFO = 3
      } else if ( LDA.LT.( K+1 ) ) {
         INFO = 6
      } else if ( INCX.EQ.0 ) {
         INFO = 8
      } else if ( INCY.EQ.0 ) {
         INFO = 11
      }
      if ( INFO.NE.0 ) {
         xerbla('ZSBMV ', INFO );
         RETURN
      }

      // Quick return if possible.

      IF( ( N.EQ.0 ) .OR. ( ( ALPHA.EQ.ZERO ) .AND. ( BETA.EQ.ONE ) ) ) RETURN

      // Set up the start points in  X  and  Y.

      if ( INCX.GT.0 ) {
         KX = 1
      } else {
         KX = 1 - ( N-1 )*INCX
      }
      if ( INCY.GT.0 ) {
         KY = 1
      } else {
         KY = 1 - ( N-1 )*INCY
      }

      // Start the operations. In this version the elements of the array A
      // are accessed sequentially with one pass through A.

      // First form  y := beta*y.

      if ( BETA.NE.ONE ) {
         if ( INCY.EQ.1 ) {
            if ( BETA.EQ.ZERO ) {
               DO 10 I = 1, N
                  Y( I ) = ZERO
   10          CONTINUE
            } else {
               DO 20 I = 1, N
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            }
         } else {
            IY = KY
            if ( BETA.EQ.ZERO ) {
               DO 30 I = 1, N
                  Y( IY ) = ZERO
                  IY = IY + INCY
   30          CONTINUE
            } else {
               DO 40 I = 1, N
                  Y( IY ) = BETA*Y( IY )
                  IY = IY + INCY
   40          CONTINUE
            }
         }
      }
      IF( ALPHA.EQ.ZERO ) RETURN
      if ( LSAME( UPLO, 'U' ) ) {

         // Form  y  when upper triangle of A is stored.

         KPLUS1 = K + 1
         if ( ( INCX.EQ.1 ) .AND. ( INCY.EQ.1 ) ) {
            DO 60 J = 1, N
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               L = KPLUS1 - J
               DO 50 I = MAX( 1, J-K ), J - 1
                  Y( I ) = Y( I ) + TEMP1*A( L+I, J )
                  TEMP2 = TEMP2 + A( L+I, J )*X( I )
   50          CONTINUE
               Y( J ) = Y( J ) + TEMP1*A( KPLUS1, J ) + ALPHA*TEMP2
   60       CONTINUE
         } else {
            JX = KX
            JY = KY
            DO 80 J = 1, N
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ZERO
               IX = KX
               IY = KY
               L = KPLUS1 - J
               DO 70 I = MAX( 1, J-K ), J - 1
                  Y( IY ) = Y( IY ) + TEMP1*A( L+I, J )
                  TEMP2 = TEMP2 + A( L+I, J )*X( IX )
                  IX = IX + INCX
                  IY = IY + INCY
   70          CONTINUE
               Y( JY ) = Y( JY ) + TEMP1*A( KPLUS1, J ) + ALPHA*TEMP2
               JX = JX + INCX
               JY = JY + INCY
               if ( J.GT.K ) {
                  KX = KX + INCX
                  KY = KY + INCY
               }
   80       CONTINUE
         }
      } else {

         // Form  y  when lower triangle of A is stored.

         if ( ( INCX.EQ.1 ) .AND. ( INCY.EQ.1 ) ) {
            DO 100 J = 1, N
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               Y( J ) = Y( J ) + TEMP1*A( 1, J )
               L = 1 - J
               DO 90 I = J + 1, MIN( N, J+K )
                  Y( I ) = Y( I ) + TEMP1*A( L+I, J )
                  TEMP2 = TEMP2 + A( L+I, J )*X( I )
   90          CONTINUE
               Y( J ) = Y( J ) + ALPHA*TEMP2
  100       CONTINUE
         } else {
            JX = KX
            JY = KY
            DO 120 J = 1, N
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ZERO
               Y( JY ) = Y( JY ) + TEMP1*A( 1, J )
               L = 1 - J
               IX = JX
               IY = JY
               DO 110 I = J + 1, MIN( N, J+K )
                  IX = IX + INCX
                  IY = IY + INCY
                  Y( IY ) = Y( IY ) + TEMP1*A( L+I, J )
                  TEMP2 = TEMP2 + A( L+I, J )*X( IX )
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP2
               JX = JX + INCX
               JY = JY + INCY
  120       CONTINUE
         }
      }

      RETURN

      // End of ZSBMV

      }
