      double sqrt17(TRANS, IRESID, M, N, NRHS, A, LDA, X, LDX, B, LDB, C, WORK, LWORK ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS;
      int                IRESID, LDA, LDB, LDX, LWORK, M, N, NRHS;
      double               A( LDA, * ), B( LDB, * ), C( LDB, * ), WORK( LWORK ), X( LDX, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                INFO, ISCL, NCOLS, NROWS;
      double               ERR, NORMA, NORMB, NORMRS, SMLNUM;
      double               RWORK( 1 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SLAMCH, SLANGE;
      // EXTERNAL lsame, SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLACPY, SLASCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, REAL

      SQRT17 = ZERO;

      if ( lsame( TRANS, 'N' ) ) {
         NROWS = M;
         NCOLS = N;
      } else if ( lsame( TRANS, 'T' ) ) {
         NROWS = N;
         NCOLS = M;
      } else {
         xerbla('SQRT17', 1 );
         return;
      }

      if ( LWORK < NCOLS*NRHS ) {
         xerbla('SQRT17', 13 );
         return;
      }

      if ( M <= 0 || N <= 0 || NRHS <= 0 ) {
         return;
      }

      NORMA = SLANGE( 'One-norm', M, N, A, LDA, RWORK );
      SMLNUM = SLAMCH( 'Safe minimum' ) / SLAMCH( 'Precision' );
      ISCL = 0;

      // compute residual and scale it

      slacpy('All', NROWS, NRHS, B, LDB, C, LDB );
      sgemm(TRANS, 'No transpose', NROWS, NRHS, NCOLS, -ONE, A, LDA, X, LDX, ONE, C, LDB );
      NORMRS = SLANGE( 'Max', NROWS, NRHS, C, LDB, RWORK );
      if ( NORMRS > SMLNUM ) {
         ISCL = 1;
         slascl('General', 0, 0, NORMRS, ONE, NROWS, NRHS, C, LDB, INFO );
      }

      // compute R**T * op(A)

      sgemm('Transpose', TRANS, NRHS, NCOLS, NROWS, ONE, C, LDB, A, LDA, ZERO, WORK, NRHS );

      // compute and properly scale error

      ERR = SLANGE( 'One-norm', NRHS, NCOLS, WORK, NRHS, RWORK );
      if (NORMA != ZERO) ERR = ERR / NORMA;

      if (ISCL == 1) ERR = ERR*NORMRS;

      if ( IRESID == 1 ) {
         NORMB = SLANGE( 'One-norm', NROWS, NRHS, B, LDB, RWORK );
         if (NORMB != ZERO) ERR = ERR / NORMB;
      } else {
         if (NORMRS != ZERO) ERR = ERR / NORMRS;
      }

      SQRT17 = ERR / ( SLAMCH( 'Epsilon' )*REAL( max( M, N, NRHS ) ) );
      }
