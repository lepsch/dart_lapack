      void dsytri2(UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             A( LDA, * ), WORK( * );
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
      // EXTERNAL DSYTRI, DSYTRI2X, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      LQUERY = ( LWORK == -1 );

      // Get blocksize

      NBMAX = ILAENV( 1, 'DSYTRI2', UPLO, N, -1, -1, -1 );
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
         xerbla('DSYTRI2', -INFO );
         return;
      } else if ( LQUERY ) {
         WORK[1] = MINSIZE;
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( NBMAX >= N ) {
         dsytri(UPLO, N, A, LDA, IPIV, WORK, INFO );
      } else {
         dsytri2x(UPLO, N, A, LDA, IPIV, WORK, NBMAX, INFO );
      }

      return;
      }