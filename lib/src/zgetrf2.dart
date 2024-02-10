      RECURSIVE SUBROUTINE ZGETRF2( M, N, A, LDA, IPIV, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, M, N;
      int                IPIV( * );
      Complex         A( LDA, * );
      // ..

      Complex         ONE, ZERO;
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ) ;
      double             SFMIN;
      Complex         TEMP;
      int                I, IINFO, N1, N2;
      // ..
      // .. External Functions ..
      //- double             DLAMCH;
      //- int                IZAMAX;
      // EXTERNAL DLAMCH, IZAMAX
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZSCAL, ZLASWP, ZTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input parameters

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('ZGETRF2', -INFO );
         return;
      }

      // Quick return if possible

      if (M == 0 || N == 0) return;

      if ( M == 1 ) {

         // Use unblocked code for one row case
         // Just need to handle IPIV and INFO

         IPIV[1] = 1;
         if ( A(1,1) == ZERO ) INFO = 1;

      } else if ( N == 1 ) {

         // Use unblocked code for one column case


         // Compute machine safe minimum

         SFMIN = dlamch('S');

         // Find pivot and test for singularity

         I = IZAMAX( M, A( 1, 1 ), 1 );
         IPIV[1] = I;
         if ( A( I, 1 ) != ZERO ) {

            // Apply the interchange

            if ( I != 1 ) {
               TEMP = A( 1, 1 );
               A[1][1] = A( I, 1 );
               A[I][1] = TEMP;
            }

            // Compute elements 2:M of the column

            if ( (A( 1, 1 )).abs() >= SFMIN ) {
               zscal(M-1, ONE / A( 1, 1 ), A( 2, 1 ), 1 );
            } else {
               for (I = 1; I <= M-1; I++) { // 10
                  A[1+I][1] = A( 1+I, 1 ) / A( 1, 1 );
               } // 10
            }

         } else {
            INFO = 1;
         }

      } else {

         // Use recursive code

         N1 = min( M, N ) / 2;
         N2 = N-N1;

                // [ A11 ]
         // Factor [ --- ]
         //        [ A21 ]

         zgetrf2(M, N1, A, LDA, IPIV, IINFO );
          if (INFO == 0 && IINFO > 0) INFO = IINFO;

                               // [ A12 ]
         // Apply interchanges to [ --- ]
         //                       [ A22 ]

         zlaswp(N2, A( 1, N1+1 ), LDA, 1, N1, IPIV, 1 );

         // Solve A12

         ztrsm('L', 'L', 'N', 'U', N1, N2, ONE, A, LDA, A( 1, N1+1 ), LDA );

         // Update A22

         zgemm('N', 'N', M-N1, N2, N1, -ONE, A( N1+1, 1 ), LDA, A( 1, N1+1 ), LDA, ONE, A( N1+1, N1+1 ), LDA );

         // Factor A22

         zgetrf2(M-N1, N2, A( N1+1, N1+1 ), LDA, IPIV( N1+1 ), IINFO );

         // Adjust INFO and the pivot indices

         if (INFO == 0 && IINFO > 0) INFO = IINFO + N1;
         for (I = N1+1; I <= min( M, N ); I++) { // 20
            IPIV[I] = IPIV( I ) + N1;
         } // 20

         // Apply interchanges to A21

         zlaswp(N1, A( 1, 1 ), LDA, N1+1, min( M, N), IPIV, 1 );

      }
      }
