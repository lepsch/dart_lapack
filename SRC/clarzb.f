      SUBROUTINE CLARZB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, V, LDV, T, LDT, C, LDC, WORK, LDWORK )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIRECT, SIDE, STOREV, TRANS;
      int                K, L, LDC, LDT, LDV, LDWORK, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            C( LDC, * ), T( LDT, * ), V( LDV, * ), WORK( LDWORK, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ONE
      const              ONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      String             TRANST;
      int                I, INFO, J;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CGEMM, CLACGV, CTRMM, XERBLA
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if (M <= 0 || N <= 0) RETURN;

      // Check for currently supported options

      INFO = 0
      if ( .NOT.LSAME( DIRECT, 'B' ) ) {
         INFO = -3
      } else if ( .NOT.LSAME( STOREV, 'R' ) ) {
         INFO = -4
      }
      if ( INFO != 0 ) {
         xerbla('CLARZB', -INFO );
         RETURN
      }

      if ( LSAME( TRANS, 'N' ) ) {
         TRANST = 'C'
      } else {
         TRANST = 'N'
      }

      if ( LSAME( SIDE, 'L' ) ) {

         // Form  H * C  or  H**H * C

         // W( 1:n, 1:k ) = C( 1:k, 1:n )**H

         for (J = 1; J <= K; J++) { // 10
            ccopy(N, C( J, 1 ), LDC, WORK( 1, J ), 1 );
         } // 10

         // W( 1:n, 1:k ) = W( 1:n, 1:k ) + ...
                         // C( m-l+1:m, 1:n )**H * V( 1:k, 1:l )**T

         if (L > 0) CALL CGEMM( 'Transpose', 'Conjugate transpose', N, K, L, ONE, C( M-L+1, 1 ), LDC, V, LDV, ONE, WORK, LDWORK );

         // W( 1:n, 1:k ) = W( 1:n, 1:k ) * T**T  or  W( 1:m, 1:k ) * T

         ctrmm('Right', 'Lower', TRANST, 'Non-unit', N, K, ONE, T, LDT, WORK, LDWORK );

         // C( 1:k, 1:n ) = C( 1:k, 1:n ) - W( 1:n, 1:k )**H

         for (J = 1; J <= N; J++) { // 30
            for (I = 1; I <= K; I++) { // 20
               C( I, J ) = C( I, J ) - WORK( J, I )
            } // 20
         } // 30

         // C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
                             // V( 1:k, 1:l )**H * W( 1:n, 1:k )**H

         if (L > 0) CALL CGEMM( 'Transpose', 'Transpose', L, N, K, -ONE, V, LDV, WORK, LDWORK, ONE, C( M-L+1, 1 ), LDC );

      } else if ( LSAME( SIDE, 'R' ) ) {

         // Form  C * H  or  C * H**H

         // W( 1:m, 1:k ) = C( 1:m, 1:k )

         for (J = 1; J <= K; J++) { // 40
            ccopy(M, C( 1, J ), 1, WORK( 1, J ), 1 );
         } // 40

         // W( 1:m, 1:k ) = W( 1:m, 1:k ) + ...
                         // C( 1:m, n-l+1:n ) * V( 1:k, 1:l )**H

         if (L > 0) CALL CGEMM( 'No transpose', 'Transpose', M, K, L, ONE, C( 1, N-L+1 ), LDC, V, LDV, ONE, WORK, LDWORK );

         // W( 1:m, 1:k ) = W( 1:m, 1:k ) * conjg( T )  or
                         // W( 1:m, 1:k ) * T**H

         for (J = 1; J <= K; J++) { // 50
            clacgv(K-J+1, T( J, J ), 1 );
         } // 50
         ctrmm('Right', 'Lower', TRANS, 'Non-unit', M, K, ONE, T, LDT, WORK, LDWORK );
         for (J = 1; J <= K; J++) { // 60
            clacgv(K-J+1, T( J, J ), 1 );
         } // 60

         // C( 1:m, 1:k ) = C( 1:m, 1:k ) - W( 1:m, 1:k )

         for (J = 1; J <= K; J++) { // 80
            for (I = 1; I <= M; I++) { // 70
               C( I, J ) = C( I, J ) - WORK( I, J )
            } // 70
         } // 80

         // C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
                             // W( 1:m, 1:k ) * conjg( V( 1:k, 1:l ) )

         for (J = 1; J <= L; J++) { // 90
            clacgv(K, V( 1, J ), 1 );
         } // 90
         if (L > 0) CALL CGEMM( 'No transpose', 'No transpose', M, L, K, -ONE, WORK, LDWORK, V, LDV, ONE, C( 1, N-L+1 ), LDC );
         for (J = 1; J <= L; J++) { // 100
            clacgv(K, V( 1, J ), 1 );
         } // 100

      }

      RETURN

      // End of CLARZB

      }
