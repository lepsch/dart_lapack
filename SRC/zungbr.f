      SUBROUTINE ZUNGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             VECT;
      int                INFO, K, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ZERO, ONE
      const              ZERO = ( 0.0D+0, 0.0D+0 ), ONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, WANTQ;
      int                I, IINFO, J, LWKOPT, MN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZUNGLQ, ZUNGQR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      WANTQ = LSAME( VECT, 'Q' )
      MN = MIN( M, N )
      LQUERY = ( LWORK == -1 )
      if ( .NOT.WANTQ .AND. .NOT.LSAME( VECT, 'P' ) ) {
         INFO = -1
      } else if ( M.LT.0 ) {
         INFO = -2
      } else if ( N.LT.0 .OR. ( WANTQ .AND. ( N.GT.M .OR. N.LT.MIN( M, K ) ) ) .OR. ( .NOT.WANTQ .AND. ( M.GT.N .OR. M.LT. MIN( N, K ) ) ) ) {
         INFO = -3
      } else if ( K.LT.0 ) {
         INFO = -4
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -6
      } else if ( LWORK.LT.MAX( 1, MN ) .AND. .NOT.LQUERY ) {
         INFO = -9
      }

      if ( INFO == 0 ) {
         WORK( 1 ) = 1
         if ( WANTQ ) {
            if ( M.GE.K ) {
               zungqr(M, N, K, A, LDA, TAU, WORK, -1, IINFO );
            } else {
               if ( M.GT.1 ) {
                  zungqr(M-1, M-1, M-1, A, LDA, TAU, WORK, -1, IINFO );
               }
            }
         } else {
            if ( K.LT.N ) {
               zunglq(M, N, K, A, LDA, TAU, WORK, -1, IINFO );
            } else {
               if ( N.GT.1 ) {
                  zunglq(N-1, N-1, N-1, A, LDA, TAU, WORK, -1, IINFO );
               }
            }
         }
         LWKOPT = INT( DBLE( WORK( 1 ) ) )
         LWKOPT = MAX (LWKOPT, MN)
      }

      if ( INFO.NE.0 ) {
         xerbla('ZUNGBR', -INFO );
         RETURN
      } else if ( LQUERY ) {
         WORK( 1 ) = LWKOPT
         RETURN
      }

      // Quick return if possible

      if ( M == 0 .OR. N == 0 ) {
         WORK( 1 ) = 1
         RETURN
      }

      if ( WANTQ ) {

         // Form Q, determined by a call to ZGEBRD to reduce an m-by-k
         // matrix

         if ( M.GE.K ) {

            // If m >= k, assume m >= n >= k

            zungqr(M, N, K, A, LDA, TAU, WORK, LWORK, IINFO );

         } else {

            // If m < k, assume m = n

            // Shift the vectors which define the elementary reflectors one
            // column to the right, and set the first row and column of Q
            // to those of the unit matrix

            DO 20 J = M, 2, -1
               A( 1, J ) = ZERO
               for (I = J + 1; I <= M; I++) { // 10
                  A( I, J ) = A( I, J-1 )
               } // 10
            } // 20
            A( 1, 1 ) = ONE
            for (I = 2; I <= M; I++) { // 30
               A( I, 1 ) = ZERO
            } // 30
            if ( M.GT.1 ) {

               // Form Q(2:m,2:m)

               zungqr(M-1, M-1, M-1, A( 2, 2 ), LDA, TAU, WORK, LWORK, IINFO );
            }
         }
      } else {

         // Form P**H, determined by a call to ZGEBRD to reduce a k-by-n
         // matrix

         if ( K.LT.N ) {

            // If k < n, assume k <= m <= n

            zunglq(M, N, K, A, LDA, TAU, WORK, LWORK, IINFO );

         } else {

            // If k >= n, assume m = n

            // Shift the vectors which define the elementary reflectors one
            // row downward, and set the first row and column of P**H to
            // those of the unit matrix

            A( 1, 1 ) = ONE
            for (I = 2; I <= N; I++) { // 40
               A( I, 1 ) = ZERO
            } // 40
            for (J = 2; J <= N; J++) { // 60
               DO 50 I = J - 1, 2, -1
                  A( I, J ) = A( I-1, J )
               } // 50
               A( 1, J ) = ZERO
            } // 60
            if ( N.GT.1 ) {

               // Form P**H(2:n,2:n)

               zunglq(N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK, LWORK, IINFO );
            }
         }
      }
      WORK( 1 ) = LWKOPT
      RETURN

      // End of ZUNGBR

      }
