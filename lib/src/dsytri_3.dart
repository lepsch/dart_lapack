      void dsytri_3(UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             A( LDA, * ), E( * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               UPPER, LQUERY;
      int                LWKOPT, NB;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      // EXTERNAL lsame, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSYTRI_3X, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      LQUERY = ( LWORK == -1 );

      // Determine the block size

      if ( N == 0 ) {
         LWKOPT = 1;
      } else {
         NB = max( 1, ILAENV( 1, 'DSYTRI_3', UPLO, N, -1, -1, -1 ) );
         LWKOPT = ( N+NB+1 ) * ( NB+3 );
      }
      WORK[1] = LWKOPT;

      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      } else if ( LWORK < LWKOPT && !LQUERY ) {
         INFO = -8;
      }

      if ( INFO != 0 ) {
         xerbla('DSYTRI_3', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      dsytri_3x(UPLO, N, A, LDA, E, IPIV, WORK, NB, INFO );

      WORK[1] = LWKOPT;

      return;
      }