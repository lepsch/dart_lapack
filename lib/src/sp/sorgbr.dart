      void sorgbr(final int VECT, final int M, final int N, final int K, final Matrix<double> A, final int LDA, final int TAU, final Array<double> WORK, final int LWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             VECT;
      int                INFO, K, LDA, LWORK, M, N;
      double               A( LDA, * ), TAU( * ), WORK( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               LQUERY, WANTQ;
      int                I, IINFO, J, LWKOPT, MN;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL lsame, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SORGLQ, SORGQR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input arguments

      INFO = 0;
      WANTQ = lsame( VECT, 'Q' );
      MN = min( M, N );
      LQUERY = ( LWORK == -1 );
      if ( !WANTQ && !lsame( VECT, 'P' ) ) {
         INFO = -1;
      } else if ( M < 0 ) {
         INFO = -2;
      } else if ( N < 0 || ( WANTQ && ( N > M || N < min( M, K ) ) ) || ( !WANTQ && ( M > N || M < min( N, K ) ) ) ) {
         INFO = -3;
      } else if ( K < 0 ) {
         INFO = -4;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -6;
      } else if ( LWORK < max( 1, MN ) && !LQUERY ) {
         INFO = -9;
      }

      if ( INFO == 0 ) {
         WORK[1] = 1;
         if ( WANTQ ) {
            if ( M >= K ) {
               sorgqr(M, N, K, A, LDA, TAU, WORK, -1, IINFO );
            } else {
               if ( M > 1 ) {
                  sorgqr(M-1, M-1, M-1, A, LDA, TAU, WORK, -1, IINFO );
               }
            }
         } else {
            if ( K < N ) {
               sorglq(M, N, K, A, LDA, TAU, WORK, -1, IINFO );
            } else {
               if ( N > 1 ) {
                  sorglq(N-1, N-1, N-1, A, LDA, TAU, WORK, -1, IINFO );
               }
            }
         }
         LWKOPT = INT( WORK( 1 ) );
         LWKOPT = max(LWKOPT, MN);
      }

      if ( INFO != 0 ) {
         xerbla('SORGBR', -INFO );
         return;
      } else if ( LQUERY ) {
         WORK[1] = SROUNDUP_LWORK(LWKOPT);
         return;
      }

      // Quick return if possible

      if ( M == 0 || N == 0 ) {
         WORK[1] = 1;
         return;
      }

      if ( WANTQ ) {

         // Form Q, determined by a call to SGEBRD to reduce an m-by-k
         // matrix

         if ( M >= K ) {

            // If m >= k, assume m >= n >= k

            sorgqr(M, N, K, A, LDA, TAU, WORK, LWORK, IINFO );

         } else {

            // If m < k, assume m = n

            // Shift the vectors which define the elementary reflectors one
            // column to the right, and set the first row and column of Q
            // to those of the unit matrix

            for (J = M; J >= 2; J--) { // 20
               A[1][J] = ZERO;
               for (I = J + 1; I <= M; I++) { // 10
                  A[I][J] = A( I, J-1 );
               } // 10
            } // 20
            A[1][1] = ONE;
            for (I = 2; I <= M; I++) { // 30
               A[I][1] = ZERO;
            } // 30
            if ( M > 1 ) {

               // Form Q(2:m,2:m)

               sorgqr(M-1, M-1, M-1, A( 2, 2 ), LDA, TAU, WORK, LWORK, IINFO );
            }
         }
      } else {

         // Form P**T, determined by a call to SGEBRD to reduce a k-by-n
         // matrix

         if ( K < N ) {

            // If k < n, assume k <= m <= n

            sorglq(M, N, K, A, LDA, TAU, WORK, LWORK, IINFO );

         } else {

            // If k >= n, assume m = n

            // Shift the vectors which define the elementary reflectors one
            // row downward, and set the first row and column of P**T to
            // those of the unit matrix

            A[1][1] = ONE;
            for (I = 2; I <= N; I++) { // 40
               A[I][1] = ZERO;
            } // 40
            for (J = 2; J <= N; J++) { // 60
               for (I = J - 1; I >= 2; I--) { // 50
                  A[I][J] = A( I-1, J );
               } // 50
               A[1][J] = ZERO;
            } // 60
            if ( N > 1 ) {

               // Form P**T(2:n,2:n)

               sorglq(N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK, LWORK, IINFO );
            }
         }
      }
      WORK[1] = SROUNDUP_LWORK(LWKOPT);
      }
