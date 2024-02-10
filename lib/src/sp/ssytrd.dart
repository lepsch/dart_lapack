      void ssytrd(UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, LWORK, N;
      double               A( LDA, * ), D( * ), E( * ), TAU( * ), WORK( * );
      // ..

      double               ONE;
      const              ONE = 1.0 ;
      bool               LQUERY, UPPER;
      int                I, IINFO, IWS, J, KK, LDWORK, LWKOPT, NB, NBMIN, NX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLATRD, SSYR2K, SSYTD2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL lsame, ILAENV, SROUNDUP_LWORK

      // Test the input parameters

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      LQUERY = ( LWORK == -1 );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      } else if ( LWORK < 1 && !LQUERY ) {
         INFO = -9;
      }

      if ( INFO == 0 ) {

         // Determine the block size.

         NB = ilaenv( 1, 'SSYTRD', UPLO, N, -1, -1, -1 );
         LWKOPT = N*NB;
         WORK[1] = SROUNDUP_LWORK(LWKOPT);
      }

      if ( INFO != 0 ) {
         xerbla('SSYTRD', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( N == 0 ) {
         WORK[1] = 1;
         return;
      }

      NX = N;
      IWS = 1;
      if ( NB > 1 && NB < N ) {

         // Determine when to cross over from blocked to unblocked code
         // (last block is always handled by unblocked code).

         NX = max( NB, ilaenv( 3, 'SSYTRD', UPLO, N, -1, -1, -1 ) );
         if ( NX < N ) {

            // Determine if workspace is large enough for blocked code.

            LDWORK = N;
            IWS = LDWORK*NB;
            if ( LWORK < IWS ) {

               // Not enough workspace to use optimal NB:  determine the
               // minimum value of NB, and reduce NB or force use of
               // unblocked code by setting NX = N.

               NB = max( LWORK / LDWORK, 1 );
               NBMIN = ilaenv( 2, 'SSYTRD', UPLO, N, -1, -1, -1 );
               if (NB < NBMIN) NX = N;
            }
         } else {
            NX = N;
         }
      } else {
         NB = 1;
      }

      if ( UPPER ) {

         // Reduce the upper triangle of A.
         // Columns 1:kk are handled by the unblocked method.

         KK = N - ( ( N-NX+NB-1 ) / NB )*NB;
         for (I = N - NB + 1; -NB < 0 ? I >= KK + 1 : I <= KK + 1; I += -NB) { // 20

            // Reduce columns i:i+nb-1 to tridiagonal form and form the
            // matrix W which is needed to update the unreduced part of
            // the matrix

            slatrd(UPLO, I+NB-1, NB, A, LDA, E, TAU, WORK, LDWORK );

            // Update the unreduced submatrix A(1:i-1,1:i-1), using an
            // update of the form:  A := A - V*W**T - W*V**T

            ssyr2k(UPLO, 'No transpose', I-1, NB, -ONE, A( 1, I ), LDA, WORK, LDWORK, ONE, A, LDA );

            // Copy superdiagonal elements back into A, and diagonal
            // elements into D

            for (J = I; J <= I + NB - 1; J++) { // 10
               A[J-1][J] = E( J-1 );
               D[J] = A( J, J );
            } // 10
         } // 20

         // Use unblocked code to reduce the last or only block

         ssytd2(UPLO, KK, A, LDA, D, E, TAU, IINFO );
      } else {

         // Reduce the lower triangle of A

         for (I = 1; NB < 0 ? I >= N - NX : I <= N - NX; I += NB) { // 40

            // Reduce columns i:i+nb-1 to tridiagonal form and form the
            // matrix W which is needed to update the unreduced part of
            // the matrix

            slatrd(UPLO, N-I+1, NB, A( I, I ), LDA, E( I ), TAU( I ), WORK, LDWORK );

            // Update the unreduced submatrix A(i+ib:n,i+ib:n), using
            // an update of the form:  A := A - V*W**T - W*V**T

            ssyr2k(UPLO, 'No transpose', N-I-NB+1, NB, -ONE, A( I+NB, I ), LDA, WORK( NB+1 ), LDWORK, ONE, A( I+NB, I+NB ), LDA );

            // Copy subdiagonal elements back into A, and diagonal
            // elements into D

            for (J = I; J <= I + NB - 1; J++) { // 30
               A[J+1][J] = E( J );
               D[J] = A( J, J );
            } // 30
         } // 40

         // Use unblocked code to reduce the last or only block

         ssytd2(UPLO, N-I+1, A( I, I ), LDA, D( I ), E( I ), TAU( I ), IINFO );
      }

      WORK[1] = SROUNDUP_LWORK(LWKOPT);
      }
