      void cspmv(UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INCX, INCY, N;
      COMPLEX            ALPHA, BETA;
      // ..
      // .. Array Arguments ..
      COMPLEX            AP( * ), X( * ), Y( * );
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX            ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      COMPLEX            ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, IX, IY, J, JX, JY, K, KK, KX, KY;
      COMPLEX            TEMP1, TEMP2;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      if ( !lsame( UPLO, 'U' ) && !lsame( UPLO, 'L' ) ) {
         INFO = 1;
      } else if ( N < 0 ) {
         INFO = 2;
      } else if ( INCX == 0 ) {
         INFO = 6;
      } else if ( INCY == 0 ) {
         INFO = 9;
      }
      if ( INFO != 0 ) {
         xerbla('CSPMV ', INFO );
         return;
      }

      // Quick return if possible.

      if( ( N == 0 ) || ( ( ALPHA == ZERO ) && ( BETA == ONE ) ) ) return;

      // Set up the start points in  X  and  Y.

      if ( INCX > 0 ) {
         KX = 1;
      } else {
         KX = 1 - ( N-1 )*INCX;
      }
      if ( INCY > 0 ) {
         KY = 1;
      } else {
         KY = 1 - ( N-1 )*INCY;
      }

      // Start the operations. In this version the elements of the array AP
      // are accessed sequentially with one pass through AP.

      // First form  y := beta*y.

      if ( BETA != ONE ) {
         if ( INCY == 1 ) {
            if ( BETA == ZERO ) {
               for (I = 1; I <= N; I++) { // 10
                  Y[I] = ZERO;
               } // 10
            } else {
               for (I = 1; I <= N; I++) { // 20
                  Y[I] = BETA*Y( I );
               } // 20
            }
         } else {
            IY = KY;
            if ( BETA == ZERO ) {
               for (I = 1; I <= N; I++) { // 30
                  Y[IY] = ZERO;
                  IY = IY + INCY;
               } // 30
            } else {
               for (I = 1; I <= N; I++) { // 40
                  Y[IY] = BETA*Y( IY );
                  IY = IY + INCY;
               } // 40
            }
         }
      }
      if (ALPHA == ZERO) return;
      KK = 1;
      if ( lsame( UPLO, 'U' ) ) {

         // Form  y  when AP contains the upper triangle.

         if ( ( INCX == 1 ) && ( INCY == 1 ) ) {
            for (J = 1; J <= N; J++) { // 60
               TEMP1 = ALPHA*X( J );
               TEMP2 = ZERO;
               K = KK;
               for (I = 1; I <= J - 1; I++) { // 50
                  Y[I] = Y( I ) + TEMP1*AP( K );
                  TEMP2 = TEMP2 + AP( K )*X( I );
                  K = K + 1;
               } // 50
               Y[J] = Y( J ) + TEMP1*AP( KK+J-1 ) + ALPHA*TEMP2;
               KK = KK + J;
            } // 60
         } else {
            JX = KX;
            JY = KY;
            for (J = 1; J <= N; J++) { // 80
               TEMP1 = ALPHA*X( JX );
               TEMP2 = ZERO;
               IX = KX;
               IY = KY;
               for (K = KK; K <= KK + J - 2; K++) { // 70
                  Y[IY] = Y( IY ) + TEMP1*AP( K );
                  TEMP2 = TEMP2 + AP( K )*X( IX );
                  IX = IX + INCX;
                  IY = IY + INCY;
               } // 70
               Y[JY] = Y( JY ) + TEMP1*AP( KK+J-1 ) + ALPHA*TEMP2;
               JX = JX + INCX;
               JY = JY + INCY;
               KK = KK + J;
            } // 80
         }
      } else {

         // Form  y  when AP contains the lower triangle.

         if ( ( INCX == 1 ) && ( INCY == 1 ) ) {
            for (J = 1; J <= N; J++) { // 100
               TEMP1 = ALPHA*X( J );
               TEMP2 = ZERO;
               Y[J] = Y( J ) + TEMP1*AP( KK );
               K = KK + 1;
               for (I = J + 1; I <= N; I++) { // 90
                  Y[I] = Y( I ) + TEMP1*AP( K );
                  TEMP2 = TEMP2 + AP( K )*X( I );
                  K = K + 1;
               } // 90
               Y[J] = Y( J ) + ALPHA*TEMP2;
               KK = KK + ( N-J+1 );
            } // 100
         } else {
            JX = KX;
            JY = KY;
            for (J = 1; J <= N; J++) { // 120
               TEMP1 = ALPHA*X( JX );
               TEMP2 = ZERO;
               Y[JY] = Y( JY ) + TEMP1*AP( KK );
               IX = JX;
               IY = JY;
               for (K = KK + 1; K <= KK + N - J; K++) { // 110
                  IX = IX + INCX;
                  IY = IY + INCY;
                  Y[IY] = Y( IY ) + TEMP1*AP( K );
                  TEMP2 = TEMP2 + AP( K )*X( IX );
               } // 110
               Y[JY] = Y( JY ) + ALPHA*TEMP2;
               JX = JX + INCX;
               JY = JY + INCY;
               KK = KK + ( N-J+1 );
            } // 120
         }
      }

      return;
      }
