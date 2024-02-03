      SUBROUTINE ZHETRI_3( UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         A( LDA, * ), E( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               UPPER, LQUERY;
      int                LWKOPT, NB;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      // EXTERNAL LSAME, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZHETRI_3X, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )

      // Determine the block size

      NB = MAX( 1, ILAENV( 1, 'ZHETRI_3', UPLO, N, -1, -1, -1 ) )
      LWKOPT = ( N+NB+1 ) * ( NB+3 )

      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      } else if ( LWORK .LT. LWKOPT .AND. .NOT.LQUERY ) {
         INFO = -8
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'ZHETRI_3', -INFO )
         RETURN
      } else if ( LQUERY ) {
         WORK( 1 ) = LWKOPT
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      CALL ZHETRI_3X( UPLO, N, A, LDA, E, IPIV, WORK, NB, INFO )

      WORK( 1 ) = LWKOPT

      RETURN

      // End of ZHETRI_3

      }
