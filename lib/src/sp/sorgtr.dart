      void sorgtr(UPLO, N, A, LDA, TAU, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, LWORK, N;
      double               A( LDA, * ), TAU( * ), WORK( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               LQUERY, UPPER;
      int                I, IINFO, J, LWKOPT, NB;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL ILAENV, lsame, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SORGQL, SORGQR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input arguments

      INFO = 0;
      LQUERY = ( LWORK == -1 );
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      } else if ( LWORK < max( 1, N-1 ) && !LQUERY ) {
         INFO = -7;
      }

      if ( INFO == 0 ) {
         if ( UPPER ) {
           NB = ilaenv( 1, 'SORGQL', ' ', N-1, N-1, N-1, -1 );
         } else {
           NB = ilaenv( 1, 'SORGQR', ' ', N-1, N-1, N-1, -1 );
         }
         LWKOPT = max( 1, N-1 )*NB;
         WORK[1] = SROUNDUP_LWORK(LWKOPT);
      }

      if ( INFO != 0 ) {
         xerbla('SORGTR', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( N == 0 ) {
         WORK[1] = 1;
         return;
      }

      if ( UPPER ) {

         // Q was determined by a call to SSYTRD with UPLO = 'U'

         // Shift the vectors which define the elementary reflectors one
         // column to the left, and set the last row and column of Q to
         // those of the unit matrix

         for (J = 1; J <= N - 1; J++) { // 20
            for (I = 1; I <= J - 1; I++) { // 10
               A[I][J] = A( I, J+1 );
            } // 10
            A[N][J] = ZERO;
         } // 20
         for (I = 1; I <= N - 1; I++) { // 30
            A[I][N] = ZERO;
         } // 30
         A[N][N] = ONE;

         // Generate Q(1:n-1,1:n-1)

         sorgql(N-1, N-1, N-1, A, LDA, TAU, WORK, LWORK, IINFO );

      } else {

         // Q was determined by a call to SSYTRD with UPLO = 'L'.

         // Shift the vectors which define the elementary reflectors one
         // column to the right, and set the first row and column of Q to
         // those of the unit matrix

         for (J = N; J >= 2; J--) { // 50
            A[1][J] = ZERO;
            for (I = J + 1; I <= N; I++) { // 40
               A[I][J] = A( I, J-1 );
            } // 40
         } // 50
         A[1][1] = ONE;
         for (I = 2; I <= N; I++) { // 60
            A[I][1] = ZERO;
         } // 60
         if ( N > 1 ) {

            // Generate Q(2:n,2:n)

            sorgqr(N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK, LWORK, IINFO );
         }
      }
      WORK[1] = SROUNDUP_LWORK(LWKOPT);
      return;
      }
