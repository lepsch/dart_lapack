      double dqrt17(TRANS, IRESID, M, N, NRHS, A, LDA, X, LDX, B, LDB, C, WORK, LWORK ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS;
      int                IRESID, LDA, LDB, LDX, LWORK, M, N, NRHS;
      double             A( LDA, * ), B( LDB, * ), C( LDB, * ), WORK( LWORK ), X( LDX, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                INFO, ISCL, NCOLS, NROWS;
      double             ERR, NORMA, NORMB, NORMRS, SMLNUM;
      double             RWORK( 1 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, DLANGE;
      // EXTERNAL lsame, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLACPY, DLASCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX

      DQRT17 = ZERO;

      if ( lsame( TRANS, 'N' ) ) {
         NROWS = M;
         NCOLS = N;
      } else if ( lsame( TRANS, 'T' ) ) {
         NROWS = N;
         NCOLS = M;
      } else {
         xerbla('DQRT17', 1 );
         return;
      }

      if ( LWORK < NCOLS*NRHS ) {
         xerbla('DQRT17', 13 );
         return;
      }

      if ( M <= 0 || N <= 0 || NRHS <= 0 ) {
         return;
      }

      NORMA = dlange( 'One-norm', M, N, A, LDA, RWORK );
      SMLNUM = dlamch( 'Safe minimum' ) / dlamch( 'Precision' );
      ISCL = 0;

      // compute residual and scale it

      dlacpy('All', NROWS, NRHS, B, LDB, C, LDB );
      dgemm(TRANS, 'No transpose', NROWS, NRHS, NCOLS, -ONE, A, LDA, X, LDX, ONE, C, LDB );
      NORMRS = dlange( 'Max', NROWS, NRHS, C, LDB, RWORK );
      if ( NORMRS > SMLNUM ) {
         ISCL = 1;
         dlascl('General', 0, 0, NORMRS, ONE, NROWS, NRHS, C, LDB, INFO );
      }

      // compute R**T * op(A)

      dgemm('Transpose', TRANS, NRHS, NCOLS, NROWS, ONE, C, LDB, A, LDA, ZERO, WORK, NRHS );

      // compute and properly scale error

      ERR = dlange( 'One-norm', NRHS, NCOLS, WORK, NRHS, RWORK );
      if (NORMA != ZERO) ERR = ERR / NORMA;

      if (ISCL == 1) ERR = ERR*NORMRS;

      if ( IRESID == 1 ) {
         NORMB = dlange( 'One-norm', NROWS, NRHS, B, LDB, RWORK );
         if (NORMB != ZERO) ERR = ERR / NORMB;
      } else {
         if (NORMRS != ZERO) ERR = ERR / NORMRS;
      }

      DQRT17 = ERR / ( dlamch( 'Epsilon' )*(max( M, N, NRHS )).toDouble() );
      return;
      }
