      SUBROUTINE ZUNGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             VECT;
      int                INFO, K, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX*16         ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
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

      INFO = 0;
      WANTQ = LSAME( VECT, 'Q' );
      MN = MIN( M, N );
      LQUERY = ( LWORK == -1 );
      if ( !WANTQ && !LSAME( VECT, 'P' ) ) {
         INFO = -1;
      } else if ( M < 0 ) {
         INFO = -2;
      } else if ( N < 0 || ( WANTQ && ( N > M || N < MIN( M, K ) ) ) || ( !WANTQ && ( M > N || M < MIN( N, K ) ) ) ) {
         INFO = -3;
      } else if ( K < 0 ) {
         INFO = -4;
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -6;
      } else if ( LWORK < MAX( 1, MN ) && !LQUERY ) {
         INFO = -9;
      }

      if ( INFO == 0 ) {
         WORK( 1 ) = 1;
         if ( WANTQ ) {
            if ( M >= K ) {
               zungqr(M, N, K, A, LDA, TAU, WORK, -1, IINFO );
            } else {
               if ( M > 1 ) {
                  zungqr(M-1, M-1, M-1, A, LDA, TAU, WORK, -1, IINFO );
               }
            }
         } else {
            if ( K < N ) {
               zunglq(M, N, K, A, LDA, TAU, WORK, -1, IINFO );
            } else {
               if ( N > 1 ) {
                  zunglq(N-1, N-1, N-1, A, LDA, TAU, WORK, -1, IINFO );
               }
            }
         }
         LWKOPT = INT( DBLE( WORK( 1 ) ) );
         LWKOPT = MAX (LWKOPT, MN);
      }

      if ( INFO != 0 ) {
         xerbla('ZUNGBR', -INFO );
         return;
      } else if ( LQUERY ) {
         WORK( 1 ) = LWKOPT;
         return;
      }

      // Quick return if possible

      if ( M == 0 || N == 0 ) {
         WORK( 1 ) = 1;
         return;
      }

      if ( WANTQ ) {

         // Form Q, determined by a call to ZGEBRD to reduce an m-by-k
         // matrix

         if ( M >= K ) {

            // If m >= k, assume m >= n >= k

            zungqr(M, N, K, A, LDA, TAU, WORK, LWORK, IINFO );

         } else {

            // If m < k, assume m = n

            // Shift the vectors which define the elementary reflectors one
            // column to the right, and set the first row and column of Q
            // to those of the unit matrix

            DO 20 J = M, 2, -1;
               A( 1, J ) = ZERO;
               for (I = J + 1; I <= M; I++) { // 10
                  A( I, J ) = A( I, J-1 );
               } // 10
            } // 20
            A( 1, 1 ) = ONE;
            for (I = 2; I <= M; I++) { // 30
               A( I, 1 ) = ZERO;
            } // 30
            if ( M > 1 ) {

               // Form Q(2:m,2:m)

               zungqr(M-1, M-1, M-1, A( 2, 2 ), LDA, TAU, WORK, LWORK, IINFO );
            }
         }
      } else {

         // Form P**H, determined by a call to ZGEBRD to reduce a k-by-n
         // matrix

         if ( K < N ) {

            // If k < n, assume k <= m <= n

            zunglq(M, N, K, A, LDA, TAU, WORK, LWORK, IINFO );

         } else {

            // If k >= n, assume m = n

            // Shift the vectors which define the elementary reflectors one
            // row downward, and set the first row and column of P**H to
            // those of the unit matrix

            A( 1, 1 ) = ONE;
            for (I = 2; I <= N; I++) { // 40
               A( I, 1 ) = ZERO;
            } // 40
            for (J = 2; J <= N; J++) { // 60
               DO 50 I = J - 1, 2, -1;
                  A( I, J ) = A( I-1, J );
               } // 50
               A( 1, J ) = ZERO;
            } // 60
            if ( N > 1 ) {

               // Form P**H(2:n,2:n)

               zunglq(N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK, LWORK, IINFO );
            }
         }
      }
      WORK( 1 ) = LWKOPT;
      return;

      // End of ZUNGBR

      }
