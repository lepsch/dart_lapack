      SUBROUTINE DORGTSQR_ROW( M, N, MB, NB, A, LDA, T, LDT, WORK, LWORK, INFO )
      IMPLICIT NONE

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int               INFO, LDA, LDT, LWORK, M, N, MB, NB;
      // ..
      // .. Array Arguments ..
      double            A( LDA, * ), T( LDT, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                NBLOCAL, MB2, M_PLUS_ONE, ITMP, IB_BOTTOM, LWORKOPT, NUM_ALL_ROW_BLOCKS, JB_T, IB, IMB, KB, KB_LAST, KNB, MB1;
      // ..
      // .. Local Arrays ..
      double             DUMMY( 1, 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARFB_GETT, DLASET, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      LQUERY  = LWORK == -1
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 || M.LT.N ) {
         INFO = -2
      } else if ( MB.LE.N ) {
         INFO = -3
      } else if ( NB.LT.1 ) {
         INFO = -4
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -6
      } else if ( LDT.LT.MAX( 1, MIN( NB, N ) ) ) {
         INFO = -8
      } else if ( LWORK.LT.1 && .NOT.LQUERY ) {
         INFO = -10
      }

      NBLOCAL = MIN( NB, N )

      // Determine the workspace size.

      if ( INFO == 0 ) {
         LWORKOPT = NBLOCAL * MAX( NBLOCAL, ( N - NBLOCAL ) )
      }

      // Handle error in the input parameters and handle the workspace query.

      if ( INFO != 0 ) {
         xerbla('DORGTSQR_ROW', -INFO );
         RETURN
      } else if ( LQUERY ) {
         WORK( 1 ) = DBLE( LWORKOPT )
         RETURN
      }

      // Quick return if possible

      if ( MIN( M, N ) == 0 ) {
         WORK( 1 ) = DBLE( LWORKOPT )
         RETURN
      }

      // (0) Set the upper-triangular part of the matrix A to zero and
      // its diagonal elements to one.

      dlaset('U', M, N, ZERO, ONE, A, LDA );

      // KB_LAST is the column index of the last column block reflector
      // in the matrices T and V.

      KB_LAST = ( ( N-1 ) / NBLOCAL ) * NBLOCAL + 1


      // (1) Bottom-up loop over row blocks of A, except the top row block.
      // NOTE: If MB>=M, then the loop is never executed.

      if ( MB.LT.M ) {

         // MB2 is the row blocking size for the row blocks before the
         // first top row block in the matrix A. IB is the row index for
         // the row blocks in the matrix A before the first top row block.
         // IB_BOTTOM is the row index for the last bottom row block
         // in the matrix A. JB_T is the column index of the corresponding
         // column block in the matrix T.

         // Initialize variables.

         // NUM_ALL_ROW_BLOCKS is the number of row blocks in the matrix A
         // including the first row block.

         MB2 = MB - N
         M_PLUS_ONE = M + 1
         ITMP = ( M - MB - 1 ) / MB2
         IB_BOTTOM = ITMP * MB2 + MB + 1
         NUM_ALL_ROW_BLOCKS = ITMP + 2
         JB_T = NUM_ALL_ROW_BLOCKS * N + 1

         DO IB = IB_BOTTOM, MB+1, -MB2

            // Determine the block size IMB for the current row block
            // in the matrix A.

            IMB = MIN( M_PLUS_ONE - IB, MB2 )

            // Determine the column index JB_T for the current column block
            // in the matrix T.

            JB_T = JB_T - N

            // Apply column blocks of H in the row block from right to left.

            // KB is the column index of the current column block reflector
            // in the matrices T and V.

            DO KB = KB_LAST, 1, -NBLOCAL

               // Determine the size of the current column block KNB in
               // the matrices T and V.

               KNB = MIN( NBLOCAL, N - KB + 1 )

               dlarfb_gett('I', IMB, N-KB+1, KNB, T( 1, JB_T+KB-1 ), LDT, A( KB, KB ), LDA, A( IB, KB ), LDA, WORK, KNB );

            }

         }

      }

      // (2) Top row block of A.
      // NOTE: If MB>=M, then we have only one row block of A of size M
      // and we work on the entire matrix A.

      MB1 = MIN( MB, M )

      // Apply column blocks of H in the top row block from right to left.

      // KB is the column index of the current block reflector in
      // the matrices T and V.

      DO KB = KB_LAST, 1, -NBLOCAL

         // Determine the size of the current column block KNB in
         // the matrices T and V.

         KNB = MIN( NBLOCAL, N - KB + 1 )

         if ( MB1-KB-KNB+1 == 0 ) {

            // In SLARFB_GETT parameters, when M=0, then the matrix B
            // does not exist, hence we need to pass a dummy array
            // reference DUMMY(1,1) to B with LDDUMMY=1.

            dlarfb_gett('N', 0, N-KB+1, KNB, T( 1, KB ), LDT, A( KB, KB ), LDA, DUMMY( 1, 1 ), 1, WORK, KNB );
         } else {
            dlarfb_gett('N', MB1-KB-KNB+1, N-KB+1, KNB, T( 1, KB ), LDT, A( KB, KB ), LDA, A( KB+KNB, KB), LDA, WORK, KNB );

         }

      }

      WORK( 1 ) = DBLE( LWORKOPT )
      RETURN

      // End of DORGTSQR_ROW

      }
