      void cgetri(N, A, LDA, IPIV, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LWORK, N;
      int                IPIV( * );
      Complex            A( LDA, * ), WORK( * );
      // ..

      Complex            ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      bool               LQUERY;
      int                I, IWS, J, JB, JJ, JP, LDWORK, LWKOPT, NB, NBMIN, NN;
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL ILAENV, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CGEMV, CSWAP, CTRSM, CTRTRI, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input parameters.

      INFO = 0;
      NB = ilaenv( 1, 'CGETRI', ' ', N, -1, -1, -1 );
      LWKOPT = max( 1, N*NB );
      WORK[1] = SROUNDUP_LWORK( LWKOPT );
      LQUERY = ( LWORK == -1 );
      if ( N < 0 ) {
         INFO = -1;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -3;
      } else if ( LWORK < max( 1, N ) && !LQUERY ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('CGETRI', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Form inv(U).  If INFO > 0 from CTRTRI, then U is singular,
      // and the inverse is not computed.

      ctrtri('Upper', 'Non-unit', N, A, LDA, INFO );
      if (INFO > 0) return;

      NBMIN = 2;
      LDWORK = N;
      if ( NB > 1 && NB < N ) {
         IWS = max( LDWORK*NB, 1 );
         if ( LWORK < IWS ) {
            NB = LWORK / LDWORK;
            NBMIN = max( 2, ilaenv( 2, 'CGETRI', ' ', N, -1, -1, -1 ) );
         }
      } else {
         IWS = N;
      }

      // Solve the equation inv(A)*L = inv(U) for inv(A).

      if ( NB < NBMIN || NB >= N ) {

         // Use unblocked code.

         for (J = N; J >= 1; J--) { // 20

            // Copy current column of L to WORK and replace with zeros.

            for (I = J + 1; I <= N; I++) { // 10
               WORK[I] = A( I, J );
               A[I][J] = ZERO;
            } // 10

            // Compute current column of inv(A).

            if (J < N) cgemv( 'No transpose', N, N-J, -ONE, A( 1, J+1 ), LDA, WORK( J+1 ), 1, ONE, A( 1, J ), 1 );
         } // 20
      } else {

         // Use blocked code.

         NN = ( ( N-1 ) / NB )*NB + 1;
         for (J = NN; -NB < 0 ? J >= 1 : J <= 1; J += -NB) { // 50
            JB = min( NB, N-J+1 );

            // Copy current block column of L to WORK and replace with
            // zeros.

            for (JJ = J; JJ <= J + JB - 1; JJ++) { // 40
               for (I = JJ + 1; I <= N; I++) { // 30
                  WORK[I+( JJ-J )*LDWORK] = A( I, JJ );
                  A[I][JJ] = ZERO;
               } // 30
            } // 40

            // Compute current block column of inv(A).

            if (J+JB <= N) cgemm( 'No transpose', 'No transpose', N, JB, N-J-JB+1, -ONE, A( 1, J+JB ), LDA, WORK( J+JB ), LDWORK, ONE, A( 1, J ), LDA );
            ctrsm('Right', 'Lower', 'No transpose', 'Unit', N, JB, ONE, WORK( J ), LDWORK, A( 1, J ), LDA );
         } // 50
      }

      // Apply column interchanges.

      for (J = N - 1; J >= 1; J--) { // 60
         JP = IPIV( J );
         if (JP != J) cswap( N, A( 1, J ), 1, A( 1, JP ), 1 );
      } // 60

      WORK[1] = SROUNDUP_LWORK( IWS );
      return;
      }
