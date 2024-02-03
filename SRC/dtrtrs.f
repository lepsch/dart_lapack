      SUBROUTINE DTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                INFO, LDA, LDB, N, NRHS;
*     ..
*     .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      bool               NOUNIT;
*     ..
*     .. External Functions ..
      bool               LSAME;
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT. LSAME( TRANS, 'T' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
*     Check for singularity.
*
      IF( NOUNIT ) THEN
         DO 10 INFO = 1, N
            IF( A( INFO, INFO ).EQ.ZERO ) RETURN
   10    CONTINUE
      END IF
      INFO = 0
*
*     Solve A * x = b  or  A**T * x = b.
*
      CALL DTRSM( 'Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B, LDB )
*
      RETURN
*
*     End of DTRTRS
*
      END
