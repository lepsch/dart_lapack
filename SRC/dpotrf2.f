      RECURSIVE SUBROUTINE DPOTRF2( UPLO, N, A, LDA, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                N1, N2, IINFO;
      // ..
      // .. External Functions ..
      bool               LSAME, DISNAN;
      // EXTERNAL LSAME, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSYRK, DTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      }
      if ( INFO != 0 ) {
         xerbla('DPOTRF2', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // N=1 case

      if ( N == 1 ) {

         // Test for non-positive-definiteness

         if ( A( 1, 1 ).LE.ZERO.OR.DISNAN( A( 1, 1 ) ) ) {
            INFO = 1
            RETURN
         }

         // Factor

         A( 1, 1 ) = SQRT( A( 1, 1 ) )

      // Use recursive code

      } else {
         N1 = N/2
         N2 = N-N1

         // Factor A11

         dpotrf2(UPLO, N1, A( 1, 1 ), LDA, IINFO );
         if ( IINFO != 0 ) {
            INFO = IINFO
            RETURN
         }

         // Compute the Cholesky factorization A = U**T*U

         if ( UPPER ) {

            // Update and scale A12

            dtrsm('L', 'U', 'T', 'N', N1, N2, ONE, A( 1, 1 ), LDA, A( 1, N1+1 ), LDA );

            // Update and factor A22

            dsyrk(UPLO, 'T', N2, N1, -ONE, A( 1, N1+1 ), LDA, ONE, A( N1+1, N1+1 ), LDA );
            dpotrf2(UPLO, N2, A( N1+1, N1+1 ), LDA, IINFO );
            if ( IINFO != 0 ) {
               INFO = IINFO + N1
               RETURN
            }

         // Compute the Cholesky factorization A = L*L**T

         } else {

            // Update and scale A21

            dtrsm('R', 'L', 'T', 'N', N2, N1, ONE, A( 1, 1 ), LDA, A( N1+1, 1 ), LDA );

            // Update and factor A22

            dsyrk(UPLO, 'N', N2, N1, -ONE, A( N1+1, 1 ), LDA, ONE, A( N1+1, N1+1 ), LDA );
            dpotrf2(UPLO, N2, A( N1+1, N1+1 ), LDA, IINFO );
            if ( IINFO != 0 ) {
               INFO = IINFO + N1
               RETURN
            }
         }
      }
      RETURN

      // End of DPOTRF2

      }
