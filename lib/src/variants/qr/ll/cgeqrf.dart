      void cgeqrf(M, N, A, LDA, TAU, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      Complex            A( LDA, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IB, IINFO, IWS, J, K, LWKOPT, NB, NBMIN, NX, LBWORK, NT, LLWORK;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEQR2, CLARFB, CLARFT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CEILING, MAX, MIN, REAL
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL ILAENV, SROUNDUP_LWORK
      // ..
      // .. Executable Statements ..

      INFO = 0;
      NBMIN = 2;
      NX = 0;
      IWS = N;
      K = min( M, N );
      NB = ILAENV( 1, 'CGEQRF', ' ', M, N, -1, -1 );

      if ( NB > 1 && NB < K ) {

         // Determine when to cross over from blocked to unblocked code.

         NX = max( 0, ILAENV( 3, 'CGEQRF', ' ', M, N, -1, -1 ) );
      }

      // Get NT, the size of the very last T, which is the left-over from in-between K-NX and K to K, eg.:

             // NB=3     2NB=6       K=10
             // |        |           |
       // 1--2--3--4--5--6--7--8--9--10
                   // |     \________/
                // K-NX=5      NT=4

      // So here 4 x 4 is the last T stored in the workspace

      NT = K-CEILING(double(K-NX)/REAL(NB))*NB;


      // optimal workspace = space for dlarfb + space for normal T's + space for the last T

      LLWORK = max(max((N-M)*K, (N-M)*NB), max(K*NB, NB*NB));
      LLWORK = CEILING(double(LLWORK)/REAL(NB));

      if ( K == 0 ) {

         LBWORK = 0;
         LWKOPT = 1;
         WORK[1] = LWKOPT;

      } else if ( NT > NB ) {

          LBWORK = K-NT;

          // Optimal workspace for dlarfb = max(1,N)*NT

          LWKOPT = (LBWORK+LLWORK)*NB;
          WORK[1] = SROUNDUP_LWORK(LWKOPT+NT*NT);

      } else {

          LBWORK = CEILING(double(K)/REAL(NB))*NB;
          LWKOPT = (LBWORK+LLWORK-NB)*NB;
          WORK[1] = SROUNDUP_LWORK(LWKOPT);

      }


      // Test the input arguments

      LQUERY = ( LWORK == -1 );
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -4;
      } else if ( !LQUERY ) {
         if( LWORK <= 0 || ( M > 0 && LWORK < max( 1, N ) ) ) INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('CGEQRF', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( K == 0 ) {
         return;
      }

      if ( NB > 1 && NB < K ) {

         if ( NX < K ) {

            // Determine if workspace is large enough for blocked code.

            if ( NT <= NB ) {
                IWS = (LBWORK+LLWORK-NB)*NB;
            } else {
                IWS = (LBWORK+LLWORK)*NB+NT*NT;
            }

            if ( LWORK < IWS ) {

               // Not enough workspace to use optimal NB:  reduce NB and
               // determine the minimum value of NB.

               if ( NT <= NB ) {
                    NB = LWORK / (LLWORK+(LBWORK-NB));
               } else {
                    NB = (LWORK-NT*NT)/(LBWORK+LLWORK);
               }
                NBMIN = max( 2, ILAENV( 2, 'CGEQRF', ' ', M, N, -1, -1 ) );
            }
         }
      }

      if ( NB >= NBMIN && NB < K && NX < K ) {

         // Use blocked code initially

         for (I = 1; NB < 0 ? I >= K - NX : I <= K - NX; I += NB) { // 10
            IB = min( K-I+1, NB );

            // Update the current column using old T's

            for (J = 1; NB < 0 ? J >= I - NB : J <= I - NB; J += NB) { // 20

               // Apply H' to A(J:M,I:I+IB-1) from the left

               clarfb('Left', 'Transpose', 'Forward', 'Columnwise', M-J+1, IB, NB, A( J, J ), LDA, WORK(J), LBWORK, A( J, I ), LDA, WORK(LBWORK*NB+NT*NT+1), IB);

20          CONTINUE

            // Compute the QR factorization of the current block
            // A(I:M,I:I+IB-1)

            cgeqr2(M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK(LBWORK*NB+NT*NT+1), IINFO );

            if ( I+IB <= N ) {

               // Form the triangular factor of the block reflector
               // H = H(i) H(i+1) . . . H(i+ib-1)

               clarft('Forward', 'Columnwise', M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK(I), LBWORK );

            }
         } // 10
      } else {
         I = 1;
      }

      // Use unblocked code to factor the last or only block.

      if ( I <= K ) {

         if ( I != 1 ) {

             for (J = 1; NB < 0 ? J >= I - NB : J <= I - NB; J += NB) { // 30

                 // Apply H' to A(J:M,I:K) from the left

                 clarfb('Left', 'Transpose', 'Forward', 'Columnwise', M-J+1, K-I+1, NB, A( J, J ), LDA, WORK(J), LBWORK, A( J, I ), LDA, WORK(LBWORK*NB+NT*NT+1), K-I+1);
30           CONTINUE
              cgeqr2(M-I+1, K-I+1, A( I, I ), LDA, TAU( I ), WORK(LBWORK*NB+NT*NT+1),IINFO );

         } else {

         // Use unblocked code to factor the last or only block.

         cgeqr2(M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK,IINFO );

         }
      }



      // Apply update to the column M+1:N when N > M

      if ( M < N && I != 1) {

          // Form the last triangular factor of the block reflector
          // H = H(i) H(i+1) . . . H(i+ib-1)

          if ( NT <= NB ) {
               clarft('Forward', 'Columnwise', M-I+1, K-I+1, A( I, I ), LDA, TAU( I ), WORK(I), LBWORK );
          } else {
               clarft('Forward', 'Columnwise', M-I+1, K-I+1, A( I, I ), LDA, TAU( I ), WORK(LBWORK*NB+1), NT );
          }


          // Apply H' to A(1:M,M+1:N) from the left

          for (J = 1; NB < 0 ? J >= K-NX : J <= K-NX; J += NB) { // 40

               IB = min( K-J+1, NB );
                clarfb('Left', 'Transpose', 'Forward', 'Columnwise', M-J+1, N-M, IB, A( J, J ), LDA, WORK(J), LBWORK, A( J, M+1 ), LDA, WORK(LBWORK*NB+NT*NT+1), N-M);

40       CONTINUE

         if ( NT <= NB ) {
             clarfb('Left', 'Transpose', 'Forward', 'Columnwise', M-J+1, N-M, K-J+1, A( J, J ), LDA, WORK(J), LBWORK, A( J, M+1 ), LDA, WORK(LBWORK*NB+NT*NT+1), N-M);
         } else {
             clarfb('Left', 'Transpose', 'Forward', 'Columnwise', M-J+1, N-M, K-J+1, A( J, J ), LDA, WORK(LBWORK*NB+1), NT, A( J, M+1 ), LDA, WORK(LBWORK*NB+NT*NT+1), N-M);
         }

      }

      WORK[1] = SROUNDUP_LWORK(IWS);
      return;
      }