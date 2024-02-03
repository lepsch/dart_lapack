      REAL             FUNCTION SQRT17( TRANS, IRESID, M, N, NRHS, A, LDA, X, LDX, B, LDB, C, WORK, LWORK )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                IRESID, LDA, LDB, LDX, LWORK, M, N, NRHS;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), C( LDB, * ), WORK( LWORK ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
      // ..
      // .. Local Scalars ..
      int                INFO, ISCL, NCOLS, NROWS;
      REAL               ERR, NORMA, NORMB, NORMRS, SMLNUM
      // ..
      // .. Local Arrays ..
      REAL               RWORK( 1 )
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, SLANGE
      // EXTERNAL LSAME, SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLACPY, SLASCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, REAL
      // ..
      // .. Executable Statements ..

      SQRT17 = ZERO

      IF( LSAME( TRANS, 'N' ) ) THEN
         NROWS = M
         NCOLS = N
      ELSE IF( LSAME( TRANS, 'T' ) ) THEN
         NROWS = N
         NCOLS = M
      ELSE
         CALL XERBLA( 'SQRT17', 1 )
         RETURN
      END IF

      IF( LWORK.LT.NCOLS*NRHS ) THEN
         CALL XERBLA( 'SQRT17', 13 )
         RETURN
      END IF

      IF( M.LE.0 .OR. N.LE.0 .OR. NRHS.LE.0 ) THEN
         RETURN
      END IF

      NORMA = SLANGE( 'One-norm', M, N, A, LDA, RWORK )
      SMLNUM = SLAMCH( 'Safe minimum' ) / SLAMCH( 'Precision' )
      ISCL = 0

      // compute residual and scale it

      CALL SLACPY( 'All', NROWS, NRHS, B, LDB, C, LDB )
      CALL SGEMM( TRANS, 'No transpose', NROWS, NRHS, NCOLS, -ONE, A, LDA, X, LDX, ONE, C, LDB )
      NORMRS = SLANGE( 'Max', NROWS, NRHS, C, LDB, RWORK )
      IF( NORMRS.GT.SMLNUM ) THEN
         ISCL = 1
         CALL SLASCL( 'General', 0, 0, NORMRS, ONE, NROWS, NRHS, C, LDB, INFO )
      END IF

      // compute R**T * op(A)

      CALL SGEMM( 'Transpose', TRANS, NRHS, NCOLS, NROWS, ONE, C, LDB, A, LDA, ZERO, WORK, NRHS )

      // compute and properly scale error

      ERR = SLANGE( 'One-norm', NRHS, NCOLS, WORK, NRHS, RWORK )
      IF( NORMA.NE.ZERO ) ERR = ERR / NORMA

      IF( ISCL.EQ.1 ) ERR = ERR*NORMRS

      IF( IRESID.EQ.1 ) THEN
         NORMB = SLANGE( 'One-norm', NROWS, NRHS, B, LDB, RWORK )
         IF( NORMB.NE.ZERO ) ERR = ERR / NORMB
      ELSE
         IF( NORMRS.NE.ZERO ) ERR = ERR / NORMRS
      END IF

      SQRT17 = ERR / ( SLAMCH( 'Epsilon' )*REAL( MAX( M, N, NRHS ) ) )
      RETURN

      // End of SQRT17

      }
