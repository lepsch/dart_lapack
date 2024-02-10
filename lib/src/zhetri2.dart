      void zhetri2(UPLO, N, final Matrix<double> A, final int LDA, final Array<int> IPIV, final Array<double> WORK, final int LWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, LWORK, N;
      int                IPIV( * );
      Complex         A( LDA, * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               UPPER, LQUERY;
      int                MINSIZE, NBMAX;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      // EXTERNAL lsame, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZHETRI2X, ZHETRI, XERBLA

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      LQUERY = ( LWORK == -1 );

      // Get blocksize

      NBMAX = ilaenv( 1, 'ZHETRF', UPLO, N, -1, -1, -1 );
      if ( N == 0 ) {
         MINSIZE = 1;
      } else if ( NBMAX >= N ) {
         MINSIZE = N;
      } else {
         MINSIZE = (N+NBMAX+1)*(NBMAX+3);
      }

      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      } else if ( LWORK < MINSIZE && !LQUERY ) {
         INFO = -7;
      }

      if ( INFO != 0 ) {
         xerbla('ZHETRI2', -INFO );
         return;
      } else if ( LQUERY ) {
         WORK[1] = MINSIZE;
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( NBMAX >= N ) {
         zhetri(UPLO, N, A, LDA, IPIV, WORK, INFO );
      } else {
         zhetri2x(UPLO, N, A, LDA, IPIV, WORK, NBMAX, INFO );
      }

      }
