      SUBROUTINE ZSYTRI2( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )

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
      // EXTERNAL ZSYTRI, ZSYTRI2X, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      // Get blocksize
      NBMAX = ILAENV( 1, 'ZSYTRI2', UPLO, N, -1, -1, -1 )
      if ( NBMAX .GE. N ) {
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
      } else if (LWORK .LT. MINSIZE .AND. .NOT.LQUERY ) {
         INFO = -7
      }

      // Quick return if possible


      if ( INFO.NE.0 ) {
         CALL XERBLA( 'ZSYTRI2', -INFO )
         RETURN
      } else if ( LQUERY ) {
         WORK(1)=MINSIZE
         RETURN
      }
      IF( N.EQ.0 ) RETURN

      if ( NBMAX .GE. N ) {
         CALL ZSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
      } else {
         CALL ZSYTRI2X( UPLO, N, A, LDA, IPIV, WORK, NBMAX, INFO )
      }
      RETURN

      // End of ZSYTRI2

      }
