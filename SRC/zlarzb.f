      SUBROUTINE ZLARZB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, V, LDV, T, LDT, C, LDC, WORK, LDWORK )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIRECT, SIDE, STOREV, TRANS;
      int                K, L, LDC, LDT, LDV, LDWORK, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         C( LDC, * ), T( LDT, * ), V( LDV, * ), WORK( LDWORK, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE
      const              ONE = ( 1.0D+0, 0.0D+0 ) ;
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
      // EXTERNAL XERBLA, ZCOPY, ZGEMM, ZLACGV, ZTRMM
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      IF( M.LE.0 .OR. N.LE.0 ) RETURN

      // Check for currently supported options

      INFO = 0
      if ( .NOT.LSAME( DIRECT, 'B' ) ) {
         INFO = -3
      } else if ( .NOT.LSAME( STOREV, 'R' ) ) {
         INFO = -4
      }
      if ( INFO.NE.0 ) {
         xerbla('ZLARZB', -INFO );
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

         DO 10 J = 1, K
            zcopy(N, C( J, 1 ), LDC, WORK( 1, J ), 1 );
   10    CONTINUE

         // W( 1:n, 1:k ) = W( 1:n, 1:k ) + ...
                         // C( m-l+1:m, 1:n )**H * V( 1:k, 1:l )**T

         IF( L.GT.0 ) CALL ZGEMM( 'Transpose', 'Conjugate transpose', N, K, L, ONE, C( M-L+1, 1 ), LDC, V, LDV, ONE, WORK, LDWORK )

         // W( 1:n, 1:k ) = W( 1:n, 1:k ) * T**T  or  W( 1:m, 1:k ) * T

         ztrmm('Right', 'Lower', TRANST, 'Non-unit', N, K, ONE, T, LDT, WORK, LDWORK );

         // C( 1:k, 1:n ) = C( 1:k, 1:n ) - W( 1:n, 1:k )**H

         DO 30 J = 1, N
            DO 20 I = 1, K
               C( I, J ) = C( I, J ) - WORK( J, I )
   20       CONTINUE
   30    CONTINUE

         // C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
                             // V( 1:k, 1:l )**H * W( 1:n, 1:k )**H

         IF( L.GT.0 ) CALL ZGEMM( 'Transpose', 'Transpose', L, N, K, -ONE, V, LDV, WORK, LDWORK, ONE, C( M-L+1, 1 ), LDC )

      } else if ( LSAME( SIDE, 'R' ) ) {

         // Form  C * H  or  C * H**H

         // W( 1:m, 1:k ) = C( 1:m, 1:k )

         DO 40 J = 1, K
            zcopy(M, C( 1, J ), 1, WORK( 1, J ), 1 );
   40    CONTINUE

         // W( 1:m, 1:k ) = W( 1:m, 1:k ) + ...
                         // C( 1:m, n-l+1:n ) * V( 1:k, 1:l )**H

         IF( L.GT.0 ) CALL ZGEMM( 'No transpose', 'Transpose', M, K, L, ONE, C( 1, N-L+1 ), LDC, V, LDV, ONE, WORK, LDWORK )

         // W( 1:m, 1:k ) = W( 1:m, 1:k ) * conjg( T )  or
                         // W( 1:m, 1:k ) * T**H

         DO 50 J = 1, K
            zlacgv(K-J+1, T( J, J ), 1 );
   50    CONTINUE
         ztrmm('Right', 'Lower', TRANS, 'Non-unit', M, K, ONE, T, LDT, WORK, LDWORK );
         DO 60 J = 1, K
            zlacgv(K-J+1, T( J, J ), 1 );
   60    CONTINUE

         // C( 1:m, 1:k ) = C( 1:m, 1:k ) - W( 1:m, 1:k )

         DO 80 J = 1, K
            DO 70 I = 1, M
               C( I, J ) = C( I, J ) - WORK( I, J )
   70       CONTINUE
   80    CONTINUE

         // C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
                             // W( 1:m, 1:k ) * conjg( V( 1:k, 1:l ) )

         DO 90 J = 1, L
            zlacgv(K, V( 1, J ), 1 );
   90    CONTINUE
         IF( L.GT.0 ) CALL ZGEMM( 'No transpose', 'No transpose', M, L, K, -ONE, WORK, LDWORK, V, LDV, ONE, C( 1, N-L+1 ), LDC )
         DO 100 J = 1, L
            zlacgv(K, V( 1, J ), 1 );
  100    CONTINUE

      }

      RETURN

      // End of ZLARZB

      }
