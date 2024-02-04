      RECURSIVE SUBROUTINE CPOTRF2( UPLO, N, A, LDA, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      COMPLEX            CONE;
      const              CONE = (1.0, 0.0) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                N1, N2, IINFO;
      REAL               AJJ;
      // ..
      // .. External Functions ..
      //- bool               LSAME, SISNAN;
      // EXTERNAL LSAME, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHERK, CTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, REAL, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('CPOTRF2', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // N=1 case

      if ( N == 1 ) {

         // Test for non-positive-definiteness

         AJJ = REAL( A( 1, 1 ) );
         if ( AJJ <= ZERO || SISNAN( AJJ ) ) {
            INFO = 1;
            return;
         }

         // Factor

         A[1, 1] = sqrt( AJJ );

      // Use recursive code

      } else {
         N1 = N/2;
         N2 = N-N1;

         // Factor A11

         cpotrf2(UPLO, N1, A( 1, 1 ), LDA, IINFO );
         if ( IINFO != 0 ) {
            INFO = IINFO;
            return;
         }

         // Compute the Cholesky factorization A = U**H*U

         if ( UPPER ) {

            // Update and scale A12

            ctrsm('L', 'U', 'C', 'N', N1, N2, CONE, A( 1, 1 ), LDA, A( 1, N1+1 ), LDA );

            // Update and factor A22

            cherk(UPLO, 'C', N2, N1, -ONE, A( 1, N1+1 ), LDA, ONE, A( N1+1, N1+1 ), LDA );

            cpotrf2(UPLO, N2, A( N1+1, N1+1 ), LDA, IINFO );

            if ( IINFO != 0 ) {
               INFO = IINFO + N1;
               return;
            }

         // Compute the Cholesky factorization A = L*L**H

         } else {

            // Update and scale A21

            ctrsm('R', 'L', 'C', 'N', N2, N1, CONE, A( 1, 1 ), LDA, A( N1+1, 1 ), LDA );

            // Update and factor A22

            cherk(UPLO, 'N', N2, N1, -ONE, A( N1+1, 1 ), LDA, ONE, A( N1+1, N1+1 ), LDA );

            cpotrf2(UPLO, N2, A( N1+1, N1+1 ), LDA, IINFO );

            if ( IINFO != 0 ) {
               INFO = IINFO + N1;
               return;
            }

         }
      }
      return;
      }
