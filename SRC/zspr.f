      SUBROUTINE ZSPR( UPLO, N, ALPHA, X, INCX, AP )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INCX, N;
      COMPLEX*16         ALPHA
      // ..
      // .. Array Arguments ..
      COMPLEX*16         AP( * ), X( * )
      // ..

* =====================================================================

      // .. Parameters ..
      COMPLEX*16         ZERO
      const              ZERO = ( 0.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, IX, J, JX, K, KK, KX;
      COMPLEX*16         TEMP
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      if ( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = 1
      } else if ( N.LT.0 ) {
         INFO = 2
      } else if ( INCX == 0 ) {
         INFO = 5
      }
      if ( INFO.NE.0 ) {
         xerbla('ZSPR  ', INFO );
         RETURN
      }

      // Quick return if possible.

      IF( ( N == 0 ) .OR. ( ALPHA == ZERO ) ) RETURN

      // Set the start point in X if the increment is not unity.

      if ( INCX.LE.0 ) {
         KX = 1 - ( N-1 )*INCX
      } else if ( INCX.NE.1 ) {
         KX = 1
      }

      // Start the operations. In this version the elements of the array AP
      // are accessed sequentially with one pass through AP.

      KK = 1
      if ( LSAME( UPLO, 'U' ) ) {

         // Form  A  when upper triangle is stored in AP.

         if ( INCX == 1 ) {
            for (J = 1; J <= N; J++) { // 20
               if ( X( J ).NE.ZERO ) {
                  TEMP = ALPHA*X( J )
                  K = KK
                  for (I = 1; I <= J - 1; I++) { // 10
                     AP( K ) = AP( K ) + X( I )*TEMP
                     K = K + 1
                  } // 10
                  AP( KK+J-1 ) = AP( KK+J-1 ) + X( J )*TEMP
               } else {
                  AP( KK+J-1 ) = AP( KK+J-1 )
               }
               KK = KK + J
            } // 20
         } else {
            JX = KX
            for (J = 1; J <= N; J++) { // 40
               if ( X( JX ).NE.ZERO ) {
                  TEMP = ALPHA*X( JX )
                  IX = KX
                  for (K = KK; K <= KK + J - 2; K++) { // 30
                     AP( K ) = AP( K ) + X( IX )*TEMP
                     IX = IX + INCX
                  } // 30
                  AP( KK+J-1 ) = AP( KK+J-1 ) + X( JX )*TEMP
               } else {
                  AP( KK+J-1 ) = AP( KK+J-1 )
               }
               JX = JX + INCX
               KK = KK + J
            } // 40
         }
      } else {

         // Form  A  when lower triangle is stored in AP.

         if ( INCX == 1 ) {
            for (J = 1; J <= N; J++) { // 60
               if ( X( J ).NE.ZERO ) {
                  TEMP = ALPHA*X( J )
                  AP( KK ) = AP( KK ) + TEMP*X( J )
                  K = KK + 1
                  for (I = J + 1; I <= N; I++) { // 50
                     AP( K ) = AP( K ) + X( I )*TEMP
                     K = K + 1
                  } // 50
               } else {
                  AP( KK ) = AP( KK )
               }
               KK = KK + N - J + 1
            } // 60
         } else {
            JX = KX
            for (J = 1; J <= N; J++) { // 80
               if ( X( JX ).NE.ZERO ) {
                  TEMP = ALPHA*X( JX )
                  AP( KK ) = AP( KK ) + TEMP*X( JX )
                  IX = JX
                  for (K = KK + 1; K <= KK + N - J; K++) { // 70
                     IX = IX + INCX
                     AP( K ) = AP( K ) + X( IX )*TEMP
                  } // 70
               } else {
                  AP( KK ) = AP( KK )
               }
               JX = JX + INCX
               KK = KK + N - J + 1
            } // 80
         }
      }

      RETURN

      // End of ZSPR

      }
