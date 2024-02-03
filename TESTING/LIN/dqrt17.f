      DOUBLE PRECISION FUNCTION DQRT17( TRANS, IRESID, M, N, NRHS, A, LDA, X, LDX, B, LDB, C, WORK, LWORK )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            IRESID, LDA, LDB, LDX, LWORK, M, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDB, * ), WORK( LWORK ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO, ISCL, NCOLS, NROWS
      DOUBLE PRECISION   ERR, NORMA, NORMB, NORMRS, SMLNUM
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   RWORK( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           LSAME, DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DLACPY, DLASCL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX
*     ..
*     .. Executable Statements ..
*
      DQRT17 = ZERO
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         NROWS = M
         NCOLS = N
      ELSE IF( LSAME( TRANS, 'T' ) ) THEN
         NROWS = N
         NCOLS = M
      ELSE
         CALL XERBLA( 'DQRT17', 1 )
         RETURN
      END IF
*
      IF( LWORK.LT.NCOLS*NRHS ) THEN
         CALL XERBLA( 'DQRT17', 13 )
         RETURN
      END IF
*
      IF( M.LE.0 .OR. N.LE.0 .OR. NRHS.LE.0 ) THEN
         RETURN
      END IF
*
      NORMA = DLANGE( 'One-norm', M, N, A, LDA, RWORK )
      SMLNUM = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      ISCL = 0
*
*     compute residual and scale it
*
      CALL DLACPY( 'All', NROWS, NRHS, B, LDB, C, LDB )
      CALL DGEMM( TRANS, 'No transpose', NROWS, NRHS, NCOLS, -ONE, A, LDA, X, LDX, ONE, C, LDB )
      NORMRS = DLANGE( 'Max', NROWS, NRHS, C, LDB, RWORK )
      IF( NORMRS.GT.SMLNUM ) THEN
         ISCL = 1
         CALL DLASCL( 'General', 0, 0, NORMRS, ONE, NROWS, NRHS, C, LDB, INFO )
      END IF
*
*     compute R**T * op(A)
*
      CALL DGEMM( 'Transpose', TRANS, NRHS, NCOLS, NROWS, ONE, C, LDB, A, LDA, ZERO, WORK, NRHS )
*
*     compute and properly scale error
*
      ERR = DLANGE( 'One-norm', NRHS, NCOLS, WORK, NRHS, RWORK )
      IF( NORMA.NE.ZERO ) ERR = ERR / NORMA
*
      IF( ISCL.EQ.1 ) ERR = ERR*NORMRS
*
      IF( IRESID.EQ.1 ) THEN
         NORMB = DLANGE( 'One-norm', NROWS, NRHS, B, LDB, RWORK )
         IF( NORMB.NE.ZERO ) ERR = ERR / NORMB
      ELSE
         IF( NORMRS.NE.ZERO ) ERR = ERR / NORMRS
      END IF
*
      DQRT17 = ERR / ( DLAMCH( 'Epsilon' )*DBLE( MAX( M, N, NRHS ) ) )
      RETURN
*
*     End of DQRT17
*
      END
