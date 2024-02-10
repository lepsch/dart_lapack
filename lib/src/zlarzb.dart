      void zlarzb(SIDE, TRANS, DIRECT, STOREV, M, N, K, L, final Matrix<double> V, final int LDV, final Matrix<double> T, final int LDT, final Matrix<double> C, final int LDC, WORK, LDWORK ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIRECT, SIDE, STOREV, TRANS;
      int                K, L, LDC, LDT, LDV, LDWORK, M, N;
      Complex         C( LDC, * ), T( LDT, * ), V( LDV, * ), WORK( LDWORK, * );
      // ..

      Complex         ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      String             TRANST;
      int                I, INFO, J;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZCOPY, ZGEMM, ZLACGV, ZTRMM

      // Quick return if possible

      if (M <= 0 || N <= 0) return;

      // Check for currently supported options

      INFO = 0;
      if ( !lsame( DIRECT, 'B' ) ) {
         INFO = -3;
      } else if ( !lsame( STOREV, 'R' ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('ZLARZB', -INFO );
         return;
      }

      if ( lsame( TRANS, 'N' ) ) {
         TRANST = 'C';
      } else {
         TRANST = 'N';
      }

      if ( lsame( SIDE, 'L' ) ) {

         // Form  H * C  or  H**H * C

         // W( 1:n, 1:k ) = C( 1:k, 1:n )**H

         for (J = 1; J <= K; J++) { // 10
            zcopy(N, C( J, 1 ), LDC, WORK( 1, J ), 1 );
         } // 10

         // W( 1:n, 1:k ) = W( 1:n, 1:k ) + ...
         //                 C( m-l+1:m, 1:n )**H * V( 1:k, 1:l )**T

         if (L > 0) zgemm( 'Transpose', 'Conjugate transpose', N, K, L, ONE, C( M-L+1, 1 ), LDC, V, LDV, ONE, WORK, LDWORK );

         // W( 1:n, 1:k ) = W( 1:n, 1:k ) * T**T  or  W( 1:m, 1:k ) * T

         ztrmm('Right', 'Lower', TRANST, 'Non-unit', N, K, ONE, T, LDT, WORK, LDWORK );

         // C( 1:k, 1:n ) = C( 1:k, 1:n ) - W( 1:n, 1:k )**H

         for (J = 1; J <= N; J++) { // 30
            for (I = 1; I <= K; I++) { // 20
               C[I][J] = C( I, J ) - WORK( J, I );
            } // 20
         } // 30

         // C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
         //                     V( 1:k, 1:l )**H * W( 1:n, 1:k )**H

         if (L > 0) zgemm( 'Transpose', 'Transpose', L, N, K, -ONE, V, LDV, WORK, LDWORK, ONE, C( M-L+1, 1 ), LDC );

      } else if ( lsame( SIDE, 'R' ) ) {

         // Form  C * H  or  C * H**H

         // W( 1:m, 1:k ) = C( 1:m, 1:k )

         for (J = 1; J <= K; J++) { // 40
            zcopy(M, C( 1, J ), 1, WORK( 1, J ), 1 );
         } // 40

         // W( 1:m, 1:k ) = W( 1:m, 1:k ) + ...
         //                 C( 1:m, n-l+1:n ) * V( 1:k, 1:l )**H

         if (L > 0) zgemm( 'No transpose', 'Transpose', M, K, L, ONE, C( 1, N-L+1 ), LDC, V, LDV, ONE, WORK, LDWORK );

         // W( 1:m, 1:k ) = W( 1:m, 1:k ) * conjg( T )  or
         //                 W( 1:m, 1:k ) * T**H

         for (J = 1; J <= K; J++) { // 50
            zlacgv(K-J+1, T( J, J ), 1 );
         } // 50
         ztrmm('Right', 'Lower', TRANS, 'Non-unit', M, K, ONE, T, LDT, WORK, LDWORK );
         for (J = 1; J <= K; J++) { // 60
            zlacgv(K-J+1, T( J, J ), 1 );
         } // 60

         // C( 1:m, 1:k ) = C( 1:m, 1:k ) - W( 1:m, 1:k )

         for (J = 1; J <= K; J++) { // 80
            for (I = 1; I <= M; I++) { // 70
               C[I][J] = C( I, J ) - WORK( I, J );
            } // 70
         } // 80

         // C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
         //                     W( 1:m, 1:k ) * conjg( V( 1:k, 1:l ) )

         for (J = 1; J <= L; J++) { // 90
            zlacgv(K, V( 1, J ), 1 );
         } // 90
         if (L > 0) zgemm( 'No transpose', 'No transpose', M, L, K, -ONE, WORK, LDWORK, V, LDV, ONE, C( 1, N-L+1 ), LDC );
         for (J = 1; J <= L; J++) { // 100
            zlacgv(K, V( 1, J ), 1 );
         } // 100

      }

      }
