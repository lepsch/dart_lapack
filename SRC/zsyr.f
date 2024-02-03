      SUBROUTINE ZSYR( UPLO, N, ALPHA, X, INCX, A, LDA )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INCX, LDA, N;
      COMPLEX*16         ALPHA
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), X( * )
      // ..

* =====================================================================

      // .. Parameters ..
      COMPLEX*16         ZERO
      const              ZERO = ( 0.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, IX, J, JX, KX;
      COMPLEX*16         TEMP
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
      if ( .NOT.LSAME( UPLO, 'U' ) && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = 1
      } else if ( N.LT.0 ) {
         INFO = 2
      } else if ( INCX == 0 ) {
         INFO = 5
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = 7
      }
      if ( INFO != 0 ) {
         xerbla('ZSYR  ', INFO );
         RETURN
      }

      // Quick return if possible.

      IF( ( N == 0 ) .OR. ( ALPHA == ZERO ) ) RETURN

      // Set the start point in X if the increment is not unity.

      if ( INCX.LE.0 ) {
         KX = 1 - ( N-1 )*INCX
      } else if ( INCX != 1 ) {
         KX = 1
      }

      // Start the operations. In this version the elements of A are
      // accessed sequentially with one pass through the triangular part
      // of A.

      if ( LSAME( UPLO, 'U' ) ) {

         // Form  A  when A is stored in upper triangle.

         if ( INCX == 1 ) {
            for (J = 1; J <= N; J++) { // 20
               if ( X( J ) != ZERO ) {
                  TEMP = ALPHA*X( J )
                  for (I = 1; I <= J; I++) { // 10
                     A( I, J ) = A( I, J ) + X( I )*TEMP
                  } // 10
               }
            } // 20
         } else {
            JX = KX
            for (J = 1; J <= N; J++) { // 40
               if ( X( JX ) != ZERO ) {
                  TEMP = ALPHA*X( JX )
                  IX = KX
                  for (I = 1; I <= J; I++) { // 30
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX = IX + INCX
                  } // 30
               }
               JX = JX + INCX
            } // 40
         }
      } else {

         // Form  A  when A is stored in lower triangle.

         if ( INCX == 1 ) {
            for (J = 1; J <= N; J++) { // 60
               if ( X( J ) != ZERO ) {
                  TEMP = ALPHA*X( J )
                  for (I = J; I <= N; I++) { // 50
                     A( I, J ) = A( I, J ) + X( I )*TEMP
                  } // 50
               }
            } // 60
         } else {
            JX = KX
            for (J = 1; J <= N; J++) { // 80
               if ( X( JX ) != ZERO ) {
                  TEMP = ALPHA*X( JX )
                  IX = JX
                  for (I = J; I <= N; I++) { // 70
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX = IX + INCX
                  } // 70
               }
               JX = JX + INCX
            } // 80
         }
      }

      RETURN

      // End of ZSYR

      }
