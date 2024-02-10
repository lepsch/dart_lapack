      void sorgtsqr_row(M, N, MB, NB, final Matrix<double> A, final int LDA, final Matrix<double> T, final int LDT, WORK, LWORK, Box<int> INFO ) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int               INFO, LDA, LDT, LWORK, M, N, MB, NB;
      double              A( LDA, * ), T( LDT, * ), WORK( * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               LQUERY;
      int                NBLOCAL, MB2, M_PLUS_ONE, ITMP, IB_BOTTOM, LWORKOPT, NUM_ALL_ROW_BLOCKS, JB_T, IB, IMB, KB, KB_LAST, KNB, MB1;
      double               DUMMY( 1, 1 );
      // ..
      // .. External Functions ..
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARFB_GETT, SLASET, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input parameters

      INFO = 0;
      LQUERY  = LWORK == -1;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 || M < N ) {
         INFO = -2;
      } else if ( MB <= N ) {
         INFO = -3;
      } else if ( NB < 1 ) {
         INFO = -4;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -6;
      } else if ( LDT < max( 1, min( NB, N ) ) ) {
         INFO = -8;
      } else if ( LWORK < 1 && !LQUERY ) {
         INFO = -10;
      }

      NBLOCAL = min( NB, N );

      // Determine the workspace size.

      if ( INFO == 0 ) {
         LWORKOPT = NBLOCAL * max( NBLOCAL, ( N - NBLOCAL ) );
      }

      // Handle error in the input parameters and handle the workspace query.

      if ( INFO != 0 ) {
         xerbla('SORGTSQR_ROW', -INFO );
         return;
      } else if ( LQUERY ) {
         WORK[1] = SROUNDUP_LWORK( LWORKOPT );
         return;
      }

      // Quick return if possible

      if ( min( M, N ) == 0 ) {
         WORK[1] = SROUNDUP_LWORK( LWORKOPT );
         return;
      }

      // (0) Set the upper-triangular part of the matrix A to zero and
      // its diagonal elements to one.

      slaset('U', M, N, ZERO, ONE, A, LDA );

      // KB_LAST is the column index of the last column block reflector
      // in the matrices T and V.

      KB_LAST = ( ( N-1 ) / NBLOCAL ) * NBLOCAL + 1;


      // (1) Bottom-up loop over row blocks of A, except the top row block.
      // NOTE: If MB>=M, then the loop is never executed.

      if ( MB < M ) {

         // MB2 is the row blocking size for the row blocks before the
         // first top row block in the matrix A. IB is the row index for
         // the row blocks in the matrix A before the first top row block.
         // IB_BOTTOM is the row index for the last bottom row block
         // in the matrix A. JB_T is the column index of the corresponding
         // column block in the matrix T.

         // Initialize variables.

         // NUM_ALL_ROW_BLOCKS is the number of row blocks in the matrix A
         // including the first row block.

         MB2 = MB - N;
         M_PLUS_ONE = M + 1;
         ITMP = ( M - MB - 1 ) / MB2;
         IB_BOTTOM = ITMP * MB2 + MB + 1;
         NUM_ALL_ROW_BLOCKS = ITMP + 2;
         JB_T = NUM_ALL_ROW_BLOCKS * N + 1;

         for (IB = IB_BOTTOM; -MB2 < 0 ? IB >= MB+1 : IB <= MB+1; IB += -MB2) {

            // Determine the block size IMB for the current row block
            // in the matrix A.

            IMB = min( M_PLUS_ONE - IB, MB2 );

            // Determine the column index JB_T for the current column block
            // in the matrix T.

            JB_T = JB_T - N;

            // Apply column blocks of H in the row block from right to left.

            // KB is the column index of the current column block reflector
            // in the matrices T and V.

            for (KB = KB_LAST; -NBLOCAL < 0 ? KB >= 1 : KB <= 1; KB += -NBLOCAL) {

               // Determine the size of the current column block KNB in
               // the matrices T and V.

               KNB = min( NBLOCAL, N - KB + 1 );

               slarfb_gett('I', IMB, N-KB+1, KNB, T( 1, JB_T+KB-1 ), LDT, A( KB, KB ), LDA, A( IB, KB ), LDA, WORK, KNB );

            }

         }

      }

      // (2) Top row block of A.
      // NOTE: If MB>=M, then we have only one row block of A of size M
      // and we work on the entire matrix A.

      MB1 = min( MB, M );

      // Apply column blocks of H in the top row block from right to left.

      // KB is the column index of the current block reflector in
      // the matrices T and V.

      for (KB = KB_LAST; -NBLOCAL < 0 ? KB >= 1 : KB <= 1; KB += -NBLOCAL) {

         // Determine the size of the current column block KNB in
         // the matrices T and V.

         KNB = min( NBLOCAL, N - KB + 1 );

         if ( MB1-KB-KNB+1 == 0 ) {

            // In SLARFB_GETT parameters, when M=0, then the matrix B
            // does not exist, hence we need to pass a dummy array
            // reference DUMMY(1,1) to B with LDDUMMY=1.

            slarfb_gett('N', 0, N-KB+1, KNB, T( 1, KB ), LDT, A( KB, KB ), LDA, DUMMY( 1, 1 ), 1, WORK, KNB );
         } else {
            slarfb_gett('N', MB1-KB-KNB+1, N-KB+1, KNB, T( 1, KB ), LDT, A( KB, KB ), LDA, A( KB+KNB, KB), LDA, WORK, KNB );

         }

      }

      WORK[1] = SROUNDUP_LWORK( LWORKOPT );
      }
