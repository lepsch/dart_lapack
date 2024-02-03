      SUBROUTINE DPPTRS( UPLO, N, NRHS, AP, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   AP( * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      bool               UPPER;
      int                I
*     ..
*     .. External Functions ..
      bool               LSAME;
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DTPSV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPPTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN
*
      IF( UPPER ) THEN
*
*        Solve A*X = B where A = U**T * U.
*
         DO 10 I = 1, NRHS
*
*           Solve U**T *X = B, overwriting B with X.
*
            CALL DTPSV( 'Upper', 'Transpose', 'Non-unit', N, AP, B( 1, I ), 1 )
*
*           Solve U*X = B, overwriting B with X.
*
            CALL DTPSV( 'Upper', 'No transpose', 'Non-unit', N, AP, B( 1, I ), 1 )
   10    CONTINUE
      ELSE
*
*        Solve A*X = B where A = L * L**T.
*
         DO 20 I = 1, NRHS
*
*           Solve L*Y = B, overwriting B with X.
*
            CALL DTPSV( 'Lower', 'No transpose', 'Non-unit', N, AP, B( 1, I ), 1 )
*
*           Solve L**T *X = Y, overwriting B with X.
*
            CALL DTPSV( 'Lower', 'Transpose', 'Non-unit', N, AP, B( 1, I ), 1 )
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of DPPTRS
*
      END
