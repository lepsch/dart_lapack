      void zhetri_3(final int UPLO, final int N, final Matrix<double> A, final int LDA, final int E, final Array<int> IPIV, final Array<double> WORK, final int LWORK, final Box<int> INFO,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, LWORK, N;
      int                IPIV( * );
      Complex         A( LDA, * ), E( * ), WORK( * );
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
      // EXTERNAL ZHETRI_3X, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      LQUERY = ( LWORK == -1 );

      // Determine the block size

      NB = max( 1, ilaenv( 1, 'ZHETRI_3', UPLO, N, -1, -1, -1 ) );
      LWKOPT = ( N+NB+1 ) * ( NB+3 );

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
         xerbla('ZHETRI_3', -INFO );
         return;
      } else if ( LQUERY ) {
         WORK[1] = LWKOPT;
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      zhetri_3x(UPLO, N, A, LDA, E, IPIV, WORK, NB, INFO );

      WORK[1] = LWKOPT;

      }
