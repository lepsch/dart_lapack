      void sgebrd(M, N, final Matrix<double> A, final int LDA, D, E, TAUQ, TAUP, final Array<double> WORK, final int LWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LWORK, M, N;
      double               A( LDA, * ), D( * ), E( * ), TAUP( * ), TAUQ( * ), WORK( * );
      // ..

      double               ONE;
      const              ONE = 1.0 ;
      bool               LQUERY;
      int                I, IINFO, J, LDWRKX, LDWRKY, LWKMIN, LWKOPT, MINMN, NB, NBMIN, NX, WS;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEBD2, SGEMM, SLABRD, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL ILAENV, SROUNDUP_LWORK

      // Test the input parameters

      INFO = 0;
      MINMN = min( M, N );
      if ( MINMN == 0 ) {
         LWKMIN = 1;
         LWKOPT = 1;
      } else {
         LWKMIN = max( M, N );
         NB = max( 1, ilaenv( 1, 'SGEBRD', ' ', M, N, -1, -1 ) );
         LWKOPT = ( M+N )*NB;
      }
      WORK[1] = SROUNDUP_LWORK( LWKOPT );
      LQUERY = ( LWORK == -1 );
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -4;
      } else if ( LWORK < LWKMIN && !LQUERY ) {
         INFO = -10;
      }
      if ( INFO < 0 ) {
         xerbla('SGEBRD', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( MINMN == 0 ) {
         WORK[1] = 1;
         return;
      }

      WS = max( M, N );
      LDWRKX = M;
      LDWRKY = N;

      if ( NB > 1 && NB < MINMN ) {

         // Set the crossover point NX.

         NX = max( NB, ilaenv( 3, 'SGEBRD', ' ', M, N, -1, -1 ) );

         // Determine when to switch from blocked to unblocked code.

         if ( NX < MINMN ) {
            WS = LWKOPT;
            if ( LWORK < WS ) {

               // Not enough work space for the optimal NB, consider using
               // a smaller block size.

               NBMIN = ilaenv( 2, 'SGEBRD', ' ', M, N, -1, -1 );
               if ( LWORK >= ( M+N )*NBMIN ) {
                  NB = LWORK / ( M+N );
               } else {
                  NB = 1;
                  NX = MINMN;
               }
            }
         }
      } else {
         NX = MINMN;
      }

      for (I = 1; NB < 0 ? I >= MINMN - NX : I <= MINMN - NX; I += NB) { // 30

         // Reduce rows and columns i:i+nb-1 to bidiagonal form and return;
         // the matrices X and Y which are needed to update the unreduced
         // part of the matrix

         slabrd(M-I+1, N-I+1, NB, A( I, I ), LDA, D( I ), E( I ), TAUQ( I ), TAUP( I ), WORK, LDWRKX, WORK( LDWRKX*NB+1 ), LDWRKY );

         // Update the trailing submatrix A(i+nb:m,i+nb:n), using an update
         // of the form  A := A - V*Y**T - X*U**T

         sgemm('No transpose', 'Transpose', M-I-NB+1, N-I-NB+1, NB, -ONE, A( I+NB, I ), LDA, WORK( LDWRKX*NB+NB+1 ), LDWRKY, ONE, A( I+NB, I+NB ), LDA );
         sgemm('No transpose', 'No transpose', M-I-NB+1, N-I-NB+1, NB, -ONE, WORK( NB+1 ), LDWRKX, A( I, I+NB ), LDA, ONE, A( I+NB, I+NB ), LDA );

         // Copy diagonal and off-diagonal elements of B back into A

         if ( M >= N ) {
            for (J = I; J <= I + NB - 1; J++) { // 10
               A[J][J] = D( J );
               A[J][J+1] = E( J );
            } // 10
         } else {
            for (J = I; J <= I + NB - 1; J++) { // 20
               A[J][J] = D( J );
               A[J+1][J] = E( J );
            } // 20
         }
      } // 30

      // Use unblocked code to reduce the remainder of the matrix

      sgebd2(M-I+1, N-I+1, A( I, I ), LDA, D( I ), E( I ), TAUQ( I ), TAUP( I ), WORK, IINFO );

      WORK[1] = SROUNDUP_LWORK( WS );
      }
