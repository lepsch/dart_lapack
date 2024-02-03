      RECURSIVE SUBROUTINE ZPOTRF2( UPLO, N, A, LDA, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      COMPLEX*16         CONE
      const              CONE = (1.0D+0, 0.0D+0) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                N1, N2, IINFO;
      double             AJJ;
      // ..
      // .. External Functions ..
      bool               LSAME, DISNAN;
      // EXTERNAL LSAME, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZHERK, ZTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, DBLE, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -4
      }
      if ( INFO != 0 ) {
         xerbla('ZPOTRF2', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // N=1 case

      if ( N == 1 ) {

         // Test for non-positive-definiteness

         AJJ = DBLE( A( 1, 1 ) )
         if ( AJJ.LE.ZERO || DISNAN( AJJ ) ) {
            INFO = 1
            RETURN
         }

         // Factor

         A( 1, 1 ) = SQRT( AJJ )

      // Use recursive code

      } else {
         N1 = N/2
         N2 = N-N1

         // Factor A11

         zpotrf2(UPLO, N1, A( 1, 1 ), LDA, IINFO );
         if ( IINFO != 0 ) {
            INFO = IINFO
            RETURN
         }

         // Compute the Cholesky factorization A = U**H*U

         if ( UPPER ) {

            // Update and scale A12

            ztrsm('L', 'U', 'C', 'N', N1, N2, CONE, A( 1, 1 ), LDA, A( 1, N1+1 ), LDA );

            // Update and factor A22

            zherk(UPLO, 'C', N2, N1, -ONE, A( 1, N1+1 ), LDA, ONE, A( N1+1, N1+1 ), LDA );
            zpotrf2(UPLO, N2, A( N1+1, N1+1 ), LDA, IINFO );
            if ( IINFO != 0 ) {
               INFO = IINFO + N1
               RETURN
            }

         // Compute the Cholesky factorization A = L*L**H

         } else {

            // Update and scale A21

            ztrsm('R', 'L', 'C', 'N', N2, N1, CONE, A( 1, 1 ), LDA, A( N1+1, 1 ), LDA );

            // Update and factor A22

            zherk(UPLO, 'N', N2, N1, -ONE, A( N1+1, 1 ), LDA, ONE, A( N1+1, N1+1 ), LDA );
            zpotrf2(UPLO, N2, A( N1+1, N1+1 ), LDA, IINFO );
            if ( IINFO != 0 ) {
               INFO = IINFO + N1
               RETURN
            }
         }
      }
      RETURN

      // End of ZPOTRF2

      }
