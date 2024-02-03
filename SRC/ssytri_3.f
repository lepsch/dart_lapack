      SUBROUTINE SSYTRI_3( UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               A( LDA, * ), E( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               UPPER, LQUERY;
      int                LWKOPT, NB;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSYTRI_3X, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK == -1 )

      // Determine the block size

      if ( N == 0 ) {
         LWKOPT = 1
      } else {
         NB = MAX( 1, ILAENV( 1, 'SSYTRI_3', UPLO, N, -1, -1, -1 ) )
         LWKOPT = ( N+NB+1 ) * ( NB+3 )
      }
      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )

      if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -4
      } else if ( LWORK < LWKOPT && .NOT.LQUERY ) {
         INFO = -8
      }

      if ( INFO != 0 ) {
         xerbla('SSYTRI_3', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      ssytri_3x(UPLO, N, A, LDA, E, IPIV, WORK, NB, INFO );

      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )

      RETURN

      // End of SSYTRI_3

      }
