      double           FUNCTION DQRT17( TRANS, IRESID, M, N, NRHS, A, LDA, X, LDX, B, LDB, C, WORK, LWORK );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                IRESID, LDA, LDB, LDX, LWORK, M, N, NRHS;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), C( LDB, * ), WORK( LWORK ), X( LDX, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      int                INFO, ISCL, NCOLS, NROWS;
      double             ERR, NORMA, NORMB, NORMRS, SMLNUM;
      // ..
      // .. Local Arrays ..
      double             RWORK( 1 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANGE;
      // EXTERNAL LSAME, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLACPY, DLASCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX
      // ..
      // .. Executable Statements ..

      DQRT17 = ZERO

      if ( LSAME( TRANS, 'N' ) ) {
         NROWS = M
         NCOLS = N
      } else if ( LSAME( TRANS, 'T' ) ) {
         NROWS = N
         NCOLS = M
      } else {
         CALL XERBLA( 'DQRT17', 1 )
         RETURN
      }

      if ( LWORK.LT.NCOLS*NRHS ) {
         CALL XERBLA( 'DQRT17', 13 )
         RETURN
      }

      if ( M.LE.0 .OR. N.LE.0 .OR. NRHS.LE.0 ) {
         RETURN
      }

      NORMA = DLANGE( 'One-norm', M, N, A, LDA, RWORK )
      SMLNUM = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      ISCL = 0

      // compute residual and scale it

      CALL DLACPY( 'All', NROWS, NRHS, B, LDB, C, LDB )
      CALL DGEMM( TRANS, 'No transpose', NROWS, NRHS, NCOLS, -ONE, A, LDA, X, LDX, ONE, C, LDB )
      NORMRS = DLANGE( 'Max', NROWS, NRHS, C, LDB, RWORK )
      if ( NORMRS.GT.SMLNUM ) {
         ISCL = 1
         CALL DLASCL( 'General', 0, 0, NORMRS, ONE, NROWS, NRHS, C, LDB, INFO )
      }

      // compute R**T * op(A)

      CALL DGEMM( 'Transpose', TRANS, NRHS, NCOLS, NROWS, ONE, C, LDB, A, LDA, ZERO, WORK, NRHS )

      // compute and properly scale error

      ERR = DLANGE( 'One-norm', NRHS, NCOLS, WORK, NRHS, RWORK )
      IF( NORMA.NE.ZERO ) ERR = ERR / NORMA

      IF( ISCL.EQ.1 ) ERR = ERR*NORMRS

      if ( IRESID.EQ.1 ) {
         NORMB = DLANGE( 'One-norm', NROWS, NRHS, B, LDB, RWORK )
         IF( NORMB.NE.ZERO ) ERR = ERR / NORMB
      } else {
         IF( NORMRS.NE.ZERO ) ERR = ERR / NORMRS
      }

      DQRT17 = ERR / ( DLAMCH( 'Epsilon' )*DBLE( MAX( M, N, NRHS ) ) )
      RETURN

      // End of DQRT17

      }
