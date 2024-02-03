      REAL             FUNCTION CQRT17( TRANS, IRESID, M, N, NRHS, A, LDA, X, LDX, B, LDB, C, WORK, LWORK );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                IRESID, LDA, LDB, LDX, LWORK, M, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), C( LDB, * ), WORK( LWORK ), X( LDX, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                INFO, ISCL, NCOLS, NROWS;
      REAL               ERR, NORMA, NORMB, NORMRS, SMLNUM;
      // ..
      // .. Local Arrays ..
      REAL               RWORK( 1 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGE, SLAMCH;
      // EXTERNAL LSAME, CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CLACPY, CLASCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, REAL
      // ..
      // .. Executable Statements ..

      CQRT17 = ZERO;

      if ( LSAME( TRANS, 'N' ) ) {
         NROWS = M;
         NCOLS = N;
      } else if ( LSAME( TRANS, 'C' ) ) {
         NROWS = N;
         NCOLS = M;
      } else {
         xerbla('CQRT17', 1 );
         RETURN;
      }

      if ( LWORK < NCOLS*NRHS ) {
         xerbla('CQRT17', 13 );
         RETURN;
      }

      if (M <= 0 || N <= 0 || NRHS <= 0) RETURN;

      NORMA = CLANGE( 'One-norm', M, N, A, LDA, RWORK );
      SMLNUM = SLAMCH( 'Safe minimum' ) / SLAMCH( 'Precision' );
      ISCL = 0;

      // compute residual and scale it

      clacpy('All', NROWS, NRHS, B, LDB, C, LDB );
      cgemm(TRANS, 'No transpose', NROWS, NRHS, NCOLS, CMPLX( -ONE ), A, LDA, X, LDX, CMPLX( ONE ), C, LDB );
      NORMRS = CLANGE( 'Max', NROWS, NRHS, C, LDB, RWORK );
      if ( NORMRS > SMLNUM ) {
         ISCL = 1;
         clascl('General', 0, 0, NORMRS, ONE, NROWS, NRHS, C, LDB, INFO );
      }

      // compute R**H * op(A)

      cgemm('Conjugate transpose', TRANS, NRHS, NCOLS, NROWS, CMPLX( ONE ), C, LDB, A, LDA, CMPLX( ZERO ), WORK, NRHS );

      // compute and properly scale error

      ERR = CLANGE( 'One-norm', NRHS, NCOLS, WORK, NRHS, RWORK );
      if (NORMA != ZERO) ERR = ERR / NORMA;

      if (ISCL == 1) ERR = ERR*NORMRS;

      if ( IRESID == 1 ) {
         NORMB = CLANGE( 'One-norm', NROWS, NRHS, B, LDB, RWORK );
         if (NORMB != ZERO) ERR = ERR / NORMB;
      } else {
         if (NORMRS != ZERO) ERR = ERR / NORMRS;
      }

      CQRT17 = ERR / ( SLAMCH( 'Epsilon' )*REAL( MAX( M, N, NRHS ) ) );
      RETURN;

      // End of CQRT17

      }
