      REAL             FUNCTION CQRT17( TRANS, IRESID, M, N, NRHS, A, LDA, X, LDX, B, LDB, C, WORK, LWORK )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             TRANS;
      int                IRESID, LDA, LDB, LDX, LWORK, M, N, NRHS;
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), C( LDB, * ), WORK( LWORK ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      int                INFO, ISCL, NCOLS, NROWS;
      REAL               ERR, NORMA, NORMB, NORMRS, SMLNUM
*     ..
*     .. Local Arrays ..
      REAL               RWORK( 1 )
*     ..
*     .. External Functions ..
      bool               LSAME;
      REAL               CLANGE, SLAMCH
      EXTERNAL           LSAME, CLANGE, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMM, CLACPY, CLASCL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, MAX, REAL
*     ..
*     .. Executable Statements ..
*
      CQRT17 = ZERO
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         NROWS = M
         NCOLS = N
      ELSE IF( LSAME( TRANS, 'C' ) ) THEN
         NROWS = N
         NCOLS = M
      ELSE
         CALL XERBLA( 'CQRT17', 1 )
         RETURN
      END IF
*
      IF( LWORK.LT.NCOLS*NRHS ) THEN
         CALL XERBLA( 'CQRT17', 13 )
         RETURN
      END IF
*
      IF( M.LE.0 .OR. N.LE.0 .OR. NRHS.LE.0 ) RETURN
*
      NORMA = CLANGE( 'One-norm', M, N, A, LDA, RWORK )
      SMLNUM = SLAMCH( 'Safe minimum' ) / SLAMCH( 'Precision' )
      ISCL = 0
*
*     compute residual and scale it
*
      CALL CLACPY( 'All', NROWS, NRHS, B, LDB, C, LDB )
      CALL CGEMM( TRANS, 'No transpose', NROWS, NRHS, NCOLS, CMPLX( -ONE ), A, LDA, X, LDX, CMPLX( ONE ), C, LDB )
      NORMRS = CLANGE( 'Max', NROWS, NRHS, C, LDB, RWORK )
      IF( NORMRS.GT.SMLNUM ) THEN
         ISCL = 1
         CALL CLASCL( 'General', 0, 0, NORMRS, ONE, NROWS, NRHS, C, LDB, INFO )
      END IF
*
*     compute R**H * op(A)
*
      CALL CGEMM( 'Conjugate transpose', TRANS, NRHS, NCOLS, NROWS, CMPLX( ONE ), C, LDB, A, LDA, CMPLX( ZERO ), WORK, NRHS )
*
*     compute and properly scale error
*
      ERR = CLANGE( 'One-norm', NRHS, NCOLS, WORK, NRHS, RWORK )
      IF( NORMA.NE.ZERO ) ERR = ERR / NORMA
*
      IF( ISCL.EQ.1 ) ERR = ERR*NORMRS
*
      IF( IRESID.EQ.1 ) THEN
         NORMB = CLANGE( 'One-norm', NROWS, NRHS, B, LDB, RWORK )
         IF( NORMB.NE.ZERO ) ERR = ERR / NORMB
      ELSE
         IF( NORMRS.NE.ZERO ) ERR = ERR / NORMRS
      END IF
*
      CQRT17 = ERR / ( SLAMCH( 'Epsilon' )*REAL( MAX( M, N, NRHS ) ) )
      RETURN
*
*     End of CQRT17
*
      END
