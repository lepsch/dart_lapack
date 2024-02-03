      RECURSIVE SUBROUTINE SLAORHR_COL_GETRFNP2( M, N, A, LDA, D, INFO );
      // IMPLICIT NONE

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), D( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      REAL               SFMIN;
      int                I, IINFO, N1, N2;
      // ..
      // .. External Functions ..
      REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SSCAL, STRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SIGN, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('SLAORHR_COL_GETRFNP2', -INFO );
         return;
      }

      // Quick return if possible

      IF( MIN( M, N ) == 0 ) RETURN;

      if ( M == 1 ) {

         // One row case, (also recursion termination case),
         // use unblocked code

         // Transfer the sign

         D( 1 ) = -SIGN( ONE, A( 1, 1 ) );

         // Construct the row of U

         A( 1, 1 ) = A( 1, 1 ) - D( 1 );

      } else if ( N == 1 ) {

         // One column case, (also recursion termination case),
         // use unblocked code

         // Transfer the sign

         D( 1 ) = -SIGN( ONE, A( 1, 1 ) );

         // Construct the row of U

         A( 1, 1 ) = A( 1, 1 ) - D( 1 );

         // Scale the elements 2:M of the column

         // Determine machine safe minimum

         SFMIN = SLAMCH('S');

         // Construct the subdiagonal elements of L

         if ( ABS( A( 1, 1 ) ) >= SFMIN ) {
            sscal(M-1, ONE / A( 1, 1 ), A( 2, 1 ), 1 );
         } else {
            for (I = 2; I <= M; I++) {
               A( I, 1 ) = A( I, 1 ) / A( 1, 1 );
            }
         }

      } else {

         // Divide the matrix B into four submatrices

         N1 = MIN( M, N ) / 2;
         N2 = N-N1;


         // Factor B11, recursive call

         slaorhr_col_getrfnp2(N1, N1, A, LDA, D, IINFO );

         // Solve for B21

         strsm('R', 'U', 'N', 'N', M-N1, N1, ONE, A, LDA, A( N1+1, 1 ), LDA );

         // Solve for B12

         strsm('L', 'L', 'N', 'U', N1, N2, ONE, A, LDA, A( 1, N1+1 ), LDA );

         // Update B22, i.e. compute the Schur complement
         // B22 := B22 - B21*B12

         sgemm('N', 'N', M-N1, N2, N1, -ONE, A( N1+1, 1 ), LDA, A( 1, N1+1 ), LDA, ONE, A( N1+1, N1+1 ), LDA );

         // Factor B22, recursive call

         slaorhr_col_getrfnp2(M-N1, N2, A( N1+1, N1+1 ), LDA, D( N1+1 ), IINFO );

      }
      return;

      // End of SLAORHR_COL_GETRFNP2

      }
