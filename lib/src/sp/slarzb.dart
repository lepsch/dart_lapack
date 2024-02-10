      void slarzb(SIDE, TRANS, DIRECT, STOREV, M, N, K, L, final Matrix<double> V, final int LDV, final Matrix<double> T, final int LDT, final Matrix<double> C, final int LDC, final Array<double> _WORK, final int LDWORK) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIRECT, SIDE, STOREV, TRANS;
      int                K, L, LDC, LDT, LDV, LDWORK, M, N;
      double               C( LDC, * ), T( LDT, * ), V( LDV, * ), WORK( LDWORK, * );
      // ..

      double               ONE;
      const              ONE = 1.0 ;
      String             TRANST;
      int                I, INFO, J;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEMM, STRMM, XERBLA

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
         xerbla('SLARZB', -INFO );
         return;
      }

      if ( lsame( TRANS, 'N' ) ) {
         TRANST = 'T';
      } else {
         TRANST = 'N';
      }

      if ( lsame( SIDE, 'L' ) ) {

         // Form  H * C  or  H**T * C

         // W( 1:n, 1:k ) = C( 1:k, 1:n )**T

         for (J = 1; J <= K; J++) { // 10
            scopy(N, C( J, 1 ), LDC, WORK( 1, J ), 1 );
         } // 10

         // W( 1:n, 1:k ) = W( 1:n, 1:k ) + ...
         //                 C( m-l+1:m, 1:n )**T * V( 1:k, 1:l )**T

         if (L > 0) sgemm( 'Transpose', 'Transpose', N, K, L, ONE, C( M-L+1, 1 ), LDC, V, LDV, ONE, WORK, LDWORK );

         // W( 1:n, 1:k ) = W( 1:n, 1:k ) * T**T  or  W( 1:m, 1:k ) * T

         strmm('Right', 'Lower', TRANST, 'Non-unit', N, K, ONE, T, LDT, WORK, LDWORK );

         // C( 1:k, 1:n ) = C( 1:k, 1:n ) - W( 1:n, 1:k )**T

         for (J = 1; J <= N; J++) { // 30
            for (I = 1; I <= K; I++) { // 20
               C[I][J] = C( I, J ) - WORK( J, I );
            } // 20
         } // 30

         // C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
         //                     V( 1:k, 1:l )**T * W( 1:n, 1:k )**T

         if (L > 0) sgemm( 'Transpose', 'Transpose', L, N, K, -ONE, V, LDV, WORK, LDWORK, ONE, C( M-L+1, 1 ), LDC );

      } else if ( lsame( SIDE, 'R' ) ) {

         // Form  C * H  or  C * H**T

         // W( 1:m, 1:k ) = C( 1:m, 1:k )

         for (J = 1; J <= K; J++) { // 40
            scopy(M, C( 1, J ), 1, WORK( 1, J ), 1 );
         } // 40

         // W( 1:m, 1:k ) = W( 1:m, 1:k ) + ...
         //                 C( 1:m, n-l+1:n ) * V( 1:k, 1:l )**T

         if (L > 0) sgemm( 'No transpose', 'Transpose', M, K, L, ONE, C( 1, N-L+1 ), LDC, V, LDV, ONE, WORK, LDWORK );

         // W( 1:m, 1:k ) = W( 1:m, 1:k ) * T  or  W( 1:m, 1:k ) * T**T

         strmm('Right', 'Lower', TRANS, 'Non-unit', M, K, ONE, T, LDT, WORK, LDWORK );

         // C( 1:m, 1:k ) = C( 1:m, 1:k ) - W( 1:m, 1:k )

         for (J = 1; J <= K; J++) { // 60
            for (I = 1; I <= M; I++) { // 50
               C[I][J] = C( I, J ) - WORK( I, J );
            } // 50
         } // 60

         // C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
         //                     W( 1:m, 1:k ) * V( 1:k, 1:l )

         if (L > 0) sgemm( 'No transpose', 'No transpose', M, L, K, -ONE, WORK, LDWORK, V, LDV, ONE, C( 1, N-L+1 ), LDC );

      }

      }
