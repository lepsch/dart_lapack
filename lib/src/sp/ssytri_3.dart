      void ssytri_3(UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, LWORK, N;
      int                IPIV( * );
      double               A( LDA, * ), E( * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               UPPER, LQUERY;
      int                LWKOPT, NB;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL lsame, ILAENV, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSYTRI_3X, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      LQUERY = ( LWORK == -1 );

      // Determine the block size

      if ( N == 0 ) {
         LWKOPT = 1;
      } else {
         NB = max( 1, ilaenv( 1, 'SSYTRI_3', UPLO, N, -1, -1, -1 ) );
         LWKOPT = ( N+NB+1 ) * ( NB+3 );
      }
      WORK[1] = SROUNDUP_LWORK( LWKOPT );

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
         xerbla('SSYTRI_3', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      ssytri_3x(UPLO, N, A, LDA, E, IPIV, WORK, NB, INFO );

      WORK[1] = SROUNDUP_LWORK( LWKOPT );

      }
