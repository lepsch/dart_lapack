      SUBROUTINE DPOTRI( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
*     ..
*     .. Array Arguments ..
      double             A( LDA, * );
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      bool               LSAME;
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAUUM, DTRTRI, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOTRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
*     Invert the triangular Cholesky factor U or L.
*
      CALL DTRTRI( UPLO, 'Non-unit', N, A, LDA, INFO )
      IF( INFO.GT.0 ) RETURN
*
*     Form inv(U) * inv(U)**T or inv(L)**T * inv(L).
*
      CALL DLAUUM( UPLO, N, A, LDA, INFO )
*
      RETURN
*
*     End of DPOTRI
*
      END
