      SUBROUTINE ZLARFB_GETT( IDENT, M, N, K, T, LDT, A, LDA, B, LDB, WORK, LDWORK );
      IMPLICIT NONE;

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             IDENT;
      int                K, LDA, LDB, LDT, LDWORK, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( LDWORK, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         CONE, CZERO;
      const              CONE = ( 1.0, 0.0 ), CZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LNOTIDENT;
      int                I, J;
      // ..
      // .. EXTERNAL FUNCTIONS ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZCOPY, ZGEMM, ZTRMM
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if (M < 0 || N <= 0 || K == 0 || K > N) RETURN;

      LNOTIDENT = !LSAME( IDENT, 'I' );

      // ------------------------------------------------------------------

      // First Step. Computation of the Column Block 2:

         // ( A2 ) := H * ( A2 )
         // ( B2 )        ( B2 )

      // ------------------------------------------------------------------

      if ( N > K ) {

         // col2_(1) Compute W2: = A2. Therefore, copy A2 = A(1:K, K+1:N)
         // into W2=WORK(1:K, 1:N-K) column-by-column.

         for (J = 1; J <= N-K; J++) {
            zcopy(K, A( 1, K+J ), 1, WORK( 1, J ), 1 );
         }

         if ( LNOTIDENT ) {

            // col2_(2) Compute W2: = (V1**H) * W2 = (A1**H) * W2,
            // V1 is not an identity matrix, but unit lower-triangular
            // V1 stored in A1 (diagonal ones are not stored).


            ztrmm('L', 'L', 'C', 'U', K, N-K, CONE, A, LDA, WORK, LDWORK );
         }

         // col2_(3) Compute W2: = W2 + (V2**H) * B2 = W2 + (B1**H) * B2
         // V2 stored in B1.

         if ( M > 0 ) {
            zgemm('C', 'N', K, N-K, M, CONE, B, LDB, B( 1, K+1 ), LDB, CONE, WORK, LDWORK );
         }

         // col2_(4) Compute W2: = T * W2,
         // T is upper-triangular.

         ztrmm('L', 'U', 'N', 'N', K, N-K, CONE, T, LDT, WORK, LDWORK );

         // col2_(5) Compute B2: = B2 - V2 * W2 = B2 - B1 * W2,
         // V2 stored in B1.

         if ( M > 0 ) {
            zgemm('N', 'N', M, N-K, K, -CONE, B, LDB, WORK, LDWORK, CONE, B( 1, K+1 ), LDB );
         }

         if ( LNOTIDENT ) {

            // col2_(6) Compute W2: = V1 * W2 = A1 * W2,
            // V1 is not an identity matrix, but unit lower-triangular,
            // V1 stored in A1 (diagonal ones are not stored).

            ztrmm('L', 'L', 'N', 'U', K, N-K, CONE, A, LDA, WORK, LDWORK );
         }

         // col2_(7) Compute A2: = A2 - W2 =
                              // = A(1:K, K+1:N-K) - WORK(1:K, 1:N-K),
         // column-by-column.

         for (J = 1; J <= N-K; J++) {
            for (I = 1; I <= K; I++) {
               A( I, K+J ) = A( I, K+J ) - WORK( I, J );
            }
         }

      }

      // ------------------------------------------------------------------

      // Second Step. Computation of the Column Block 1:

         // ( A1 ) := H * ( A1 )
         // ( B1 )        (  0 )

      // ------------------------------------------------------------------

      // col1_(1) Compute W1: = A1. Copy the upper-triangular
      // A1 = A(1:K, 1:K) into the upper-triangular
      // W1 = WORK(1:K, 1:K) column-by-column.

      for (J = 1; J <= K; J++) {
         zcopy(J, A( 1, J ), 1, WORK( 1, J ), 1 );
      }

      // Set the subdiagonal elements of W1 to zero column-by-column.

      for (J = 1; J <= K - 1; J++) {
         for (I = J + 1; I <= K; I++) {
            WORK( I, J ) = CZERO;
         }
      }

      if ( LNOTIDENT ) {

         // col1_(2) Compute W1: = (V1**H) * W1 = (A1**H) * W1,
         // V1 is not an identity matrix, but unit lower-triangular
         // V1 stored in A1 (diagonal ones are not stored),
         // W1 is upper-triangular with zeroes below the diagonal.

         ztrmm('L', 'L', 'C', 'U', K, K, CONE, A, LDA, WORK, LDWORK );
      }

      // col1_(3) Compute W1: = T * W1,
      // T is upper-triangular,
      // W1 is upper-triangular with zeroes below the diagonal.

      ztrmm('L', 'U', 'N', 'N', K, K, CONE, T, LDT, WORK, LDWORK );

      // col1_(4) Compute B1: = - V2 * W1 = - B1 * W1,
      // V2 = B1, W1 is upper-triangular with zeroes below the diagonal.

      if ( M > 0 ) {
         ztrmm('R', 'U', 'N', 'N', M, K, -CONE, WORK, LDWORK, B, LDB );
      }

      if ( LNOTIDENT ) {

         // col1_(5) Compute W1: = V1 * W1 = A1 * W1,
         // V1 is not an identity matrix, but unit lower-triangular
         // V1 stored in A1 (diagonal ones are not stored),
         // W1 is upper-triangular on input with zeroes below the diagonal,
         // and square on output.

         ztrmm('L', 'L', 'N', 'U', K, K, CONE, A, LDA, WORK, LDWORK );

         // col1_(6) Compute A1: = A1 - W1 = A(1:K, 1:K) - WORK(1:K, 1:K)
         // column-by-column. A1 is upper-triangular on input.
         // If IDENT, A1 is square on output, and W1 is square,
         // if NOT IDENT, A1 is upper-triangular on output,
         // W1 is upper-triangular.

         // col1_(6)_a Compute elements of A1 below the diagonal.

         for (J = 1; J <= K - 1; J++) {
            for (I = J + 1; I <= K; I++) {
               A( I, J ) = - WORK( I, J );
            }
         }

      }

      // col1_(6)_b Compute elements of A1 on and above the diagonal.

      for (J = 1; J <= K; J++) {
         for (I = 1; I <= J; I++) {
            A( I, J ) = A( I, J ) - WORK( I, J );
         }
      }

      RETURN;

      // End of ZLARFB_GETT

      }
