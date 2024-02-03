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
      LQUERY = ( LWORK.EQ.-1 )

      // Get blocksize

      NBMAX = ILAENV( 1, 'ZHETRF', UPLO, N, -1, -1, -1 )
      IF( N.EQ.0 ) THEN
         MINSIZE = 1
      ELSE IF( NBMAX.GE.N ) THEN
         MINSIZE = N
      ELSE
         MINSIZE = (N+NBMAX+1)*(NBMAX+3)
      END IF

      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MINSIZE .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHETRI2', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         WORK( 1 ) = MINSIZE
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      IF( NBMAX.GE.N ) THEN
         CALL ZHETRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
      ELSE
         CALL ZHETRI2X( UPLO, N, A, LDA, IPIV, WORK, NBMAX, INFO )
      END IF

      RETURN

      // End of ZHETRI2

      END
