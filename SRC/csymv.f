      SUBROUTINE CSYMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INCX, INCY, LDA, N;
      COMPLEX            ALPHA, BETA
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), X( * ), Y( * )
      // ..

* =====================================================================

      // .. Parameters ..
      COMPLEX            ONE
      const              ONE = ( 1.0E+0, 0.0E+0 ) ;
      COMPLEX            ZERO
      const              ZERO = ( 0.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, IX, IY, J, JX, JY, KX, KY;
      COMPLEX            TEMP1, TEMP2
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      if ( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) {
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
         xerbla('CSYMV ', INFO );
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

      // Start the operations. In this version the elements of A are
      // accessed sequentially with one pass through the triangular part
      // of A.

      // First form  y := beta*y.

      if ( BETA.NE.ONE ) {
         if ( INCY.EQ.1 ) {
            if ( BETA.EQ.ZERO ) {
               for (I = 1; I <= N; I++) { // 10
                  Y( I ) = ZERO
   10          CONTINUE
            } else {
               for (I = 1; I <= N; I++) { // 20
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            }
         } else {
            IY = KY
            if ( BETA.EQ.ZERO ) {
               for (I = 1; I <= N; I++) { // 30
                  Y( IY ) = ZERO
                  IY = IY + INCY
   30          CONTINUE
            } else {
               for (I = 1; I <= N; I++) { // 40
                  Y( IY ) = BETA*Y( IY )
                  IY = IY + INCY
   40          CONTINUE
            }
         }
      }
      IF( ALPHA.EQ.ZERO ) RETURN
      if ( LSAME( UPLO, 'U' ) ) {

         // Form  y  when A is stored in upper triangle.

         if ( ( INCX.EQ.1 ) .AND. ( INCY.EQ.1 ) ) {
            for (J = 1; J <= N; J++) { // 60
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               DO 50 I = 1, J - 1
                  Y( I ) = Y( I ) + TEMP1*A( I, J )
                  TEMP2 = TEMP2 + A( I, J )*X( I )
   50          CONTINUE
               Y( J ) = Y( J ) + TEMP1*A( J, J ) + ALPHA*TEMP2
   60       CONTINUE
         } else {
            JX = KX
            JY = KY
            for (J = 1; J <= N; J++) { // 80
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ZERO
               IX = KX
               IY = KY
               DO 70 I = 1, J - 1
                  Y( IY ) = Y( IY ) + TEMP1*A( I, J )
                  TEMP2 = TEMP2 + A( I, J )*X( IX )
                  IX = IX + INCX
                  IY = IY + INCY
   70          CONTINUE
               Y( JY ) = Y( JY ) + TEMP1*A( J, J ) + ALPHA*TEMP2
               JX = JX + INCX
               JY = JY + INCY
   80       CONTINUE
         }
      } else {

         // Form  y  when A is stored in lower triangle.

         if ( ( INCX.EQ.1 ) .AND. ( INCY.EQ.1 ) ) {
            for (J = 1; J <= N; J++) { // 100
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               Y( J ) = Y( J ) + TEMP1*A( J, J )
               DO 90 I = J + 1, N
                  Y( I ) = Y( I ) + TEMP1*A( I, J )
                  TEMP2 = TEMP2 + A( I, J )*X( I )
   90          CONTINUE
               Y( J ) = Y( J ) + ALPHA*TEMP2
  100       CONTINUE
         } else {
            JX = KX
            JY = KY
            for (J = 1; J <= N; J++) { // 120
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ZERO
               Y( JY ) = Y( JY ) + TEMP1*A( J, J )
               IX = JX
               IY = JY
               DO 110 I = J + 1, N
                  IX = IX + INCX
                  IY = IY + INCY
                  Y( IY ) = Y( IY ) + TEMP1*A( I, J )
                  TEMP2 = TEMP2 + A( I, J )*X( IX )
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP2
               JX = JX + INCX
               JY = JY + INCY
  120       CONTINUE
         }
      }

      RETURN

      // End of CSYMV

      }
