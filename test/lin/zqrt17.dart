      double zqrt17(TRANS, IRESID, M, N, NRHS, A, LDA, X, LDX, B, LDB, C, WORK, LWORK ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS;
      int                IRESID, LDA, LDB, LDX, LWORK, M, N, NRHS;
      Complex         A( LDA, * ), B( LDB, * ), C( LDB, * ), WORK( LWORK ), X( LDX, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                INFO, ISCL, NCOLS, NROWS;
      double             ERR, NORMA, NORMB, NORMRS, SMLNUM;
      double             RWORK( 1 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, ZLANGE;
      // EXTERNAL lsame, DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEMM, ZLACPY, ZLASCL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX

      ZQRT17 = ZERO;

      if ( lsame( TRANS, 'N' ) ) {
         NROWS = M;
         NCOLS = N;
      } else if ( lsame( TRANS, 'C' ) ) {
         NROWS = N;
         NCOLS = M;
      } else {
         xerbla('ZQRT17', 1 );
         return;
      }

      if ( LWORK < NCOLS*NRHS ) {
         xerbla('ZQRT17', 13 );
         return;
      }

      if (M <= 0 || N <= 0 || NRHS <= 0) return;

      NORMA = ZLANGE( 'One-norm', M, N, A, LDA, RWORK );
      SMLNUM = dlamch( 'Safe minimum' ) / dlamch( 'Precision' );
      ISCL = 0;

      // compute residual and scale it

      zlacpy('All', NROWS, NRHS, B, LDB, C, LDB );
      zgemm(TRANS, 'No transpose', NROWS, NRHS, NCOLS, DCMPLX( -ONE ), A, LDA, X, LDX, DCMPLX( ONE ), C, LDB );
      NORMRS = ZLANGE( 'Max', NROWS, NRHS, C, LDB, RWORK );
      if ( NORMRS > SMLNUM ) {
         ISCL = 1;
         zlascl('General', 0, 0, NORMRS, ONE, NROWS, NRHS, C, LDB, INFO );
      }

      // compute R**H * op(A)

      zgemm('Conjugate transpose', TRANS, NRHS, NCOLS, NROWS, DCMPLX( ONE ), C, LDB, A, LDA, DCMPLX( ZERO ), WORK, NRHS );

      // compute and properly scale error

      ERR = ZLANGE( 'One-norm', NRHS, NCOLS, WORK, NRHS, RWORK );
      if (NORMA != ZERO) ERR = ERR / NORMA;

      if (ISCL == 1) ERR = ERR*NORMRS;

      if ( IRESID == 1 ) {
         NORMB = ZLANGE( 'One-norm', NROWS, NRHS, B, LDB, RWORK );
         if (NORMB != ZERO) ERR = ERR / NORMB;
      } else {
         if (NORMRS != ZERO) ERR = ERR / NORMRS;
      }

      ZQRT17 = ERR / ( dlamch( 'Epsilon' )*(max( M, N, NRHS )).toDouble() );
      }
