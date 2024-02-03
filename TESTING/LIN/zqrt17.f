      double           FUNCTION ZQRT17( TRANS, IRESID, M, N, NRHS, A, LDA, X, LDX, B, LDB, C, WORK, LWORK );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                IRESID, LDA, LDB, LDX, LWORK, M, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDB, * ), WORK( LWORK ), X( LDX, * )
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
      double             DLAMCH, ZLANGE;
      // EXTERNAL LSAME, DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEMM, ZLACPY, ZLASCL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX
      // ..
      // .. Executable Statements ..

      ZQRT17 = ZERO

      if ( LSAME( TRANS, 'N' ) ) {
         NROWS = M
         NCOLS = N
      } else if ( LSAME( TRANS, 'C' ) ) {
         NROWS = N
         NCOLS = M
      } else {
         CALL XERBLA( 'ZQRT17', 1 )
         RETURN
      }

      if ( LWORK.LT.NCOLS*NRHS ) {
         CALL XERBLA( 'ZQRT17', 13 )
         RETURN
      }

      IF( M.LE.0 .OR. N.LE.0 .OR. NRHS.LE.0 ) RETURN

      NORMA = ZLANGE( 'One-norm', M, N, A, LDA, RWORK )
      SMLNUM = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      ISCL = 0

      // compute residual and scale it

      CALL ZLACPY( 'All', NROWS, NRHS, B, LDB, C, LDB )
      CALL ZGEMM( TRANS, 'No transpose', NROWS, NRHS, NCOLS, DCMPLX( -ONE ), A, LDA, X, LDX, DCMPLX( ONE ), C, LDB )
      NORMRS = ZLANGE( 'Max', NROWS, NRHS, C, LDB, RWORK )
      if ( NORMRS.GT.SMLNUM ) {
         ISCL = 1
         CALL ZLASCL( 'General', 0, 0, NORMRS, ONE, NROWS, NRHS, C, LDB, INFO )
      }

      // compute R**H * op(A)

      CALL ZGEMM( 'Conjugate transpose', TRANS, NRHS, NCOLS, NROWS, DCMPLX( ONE ), C, LDB, A, LDA, DCMPLX( ZERO ), WORK, NRHS )

      // compute and properly scale error

      ERR = ZLANGE( 'One-norm', NRHS, NCOLS, WORK, NRHS, RWORK )
      IF( NORMA.NE.ZERO ) ERR = ERR / NORMA

      IF( ISCL.EQ.1 ) ERR = ERR*NORMRS

      if ( IRESID.EQ.1 ) {
         NORMB = ZLANGE( 'One-norm', NROWS, NRHS, B, LDB, RWORK )
         IF( NORMB.NE.ZERO ) ERR = ERR / NORMB
      } else {
         IF( NORMRS.NE.ZERO ) ERR = ERR / NORMRS
      }

      ZQRT17 = ERR / ( DLAMCH( 'Epsilon' )*DBLE( MAX( M, N, NRHS ) ) )
      RETURN

      // End of ZQRT17

      }
