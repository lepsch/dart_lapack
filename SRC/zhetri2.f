      SUBROUTINE ZHETRI2( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         A( LDA, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               UPPER, LQUERY;
      int                MINSIZE, NBMAX;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      // EXTERNAL LSAME, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZHETRI2X, ZHETRI, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK == -1 )

      // Get blocksize

      NBMAX = ILAENV( 1, 'ZHETRF', UPLO, N, -1, -1, -1 )
      if ( N == 0 ) {
         MINSIZE = 1
      } else if ( NBMAX.GE.N ) {
         MINSIZE = N
      } else {
         MINSIZE = (N+NBMAX+1)*(NBMAX+3)
      }

      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      } else if ( LWORK.LT.MINSIZE .AND. .NOT.LQUERY ) {
         INFO = -7
      }

      if ( INFO.NE.0 ) {
         xerbla('ZHETRI2', -INFO );
         RETURN
      } else if ( LQUERY ) {
         WORK( 1 ) = MINSIZE
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      if ( NBMAX.GE.N ) {
         zhetri(UPLO, N, A, LDA, IPIV, WORK, INFO );
      } else {
         zhetri2x(UPLO, N, A, LDA, IPIV, WORK, NBMAX, INFO );
      }

      RETURN

      // End of ZHETRI2

      }
