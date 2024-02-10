      double cqrt17(TRANS, IRESID, M, N, NRHS, final Matrix<double> A, final int LDA, final Matrix<double> X, final int LDX, final Matrix<double> B, final int LDB, C, WORK, LWORK ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS;
      int                IRESID, LDA, LDB, LDX, LWORK, M, N, NRHS;
      Complex            A( LDA, * ), B( LDB, * ), C( LDB, * ), WORK( LWORK ), X( LDX, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                INFO, ISCL, NCOLS, NROWS;
      double               ERR, NORMA, NORMB, NORMRS, SMLNUM;
      double               RWORK( 1 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANGE, SLAMCH;
      // EXTERNAL lsame, CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CLACPY, CLASCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, REAL

      CQRT17 = ZERO;

      if ( lsame( TRANS, 'N' ) ) {
         NROWS = M;
         NCOLS = N;
      } else if ( lsame( TRANS, 'C' ) ) {
         NROWS = N;
         NCOLS = M;
      } else {
         xerbla('CQRT17', 1 );
         return;
      }

      if ( LWORK < NCOLS*NRHS ) {
         xerbla('CQRT17', 13 );
         return;
      }

      if (M <= 0 || N <= 0 || NRHS <= 0) return;

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

      CQRT17 = ERR / ( SLAMCH( 'Epsilon' )*REAL( max( M, N, NRHS ) ) );
      }
