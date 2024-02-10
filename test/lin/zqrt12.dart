      double zqrt12(M, N, final Matrix<double> A, final int LDA, S, final Array<double> WORK, final int LWORK, final Array<double> RWORK) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LWORK, M, N;
      double             RWORK( * ), S( * );
      Complex         A( LDA, * ), WORK( LWORK );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, INFO, ISCL, J, MN;
      double             ANRM, BIGNUM, NRMSVL, SMLNUM;
      double             DUMMY( 1 );
      // ..
      // .. External Functions ..
      //- double             DASUM, DLAMCH, DNRM2, ZLANGE;
      // EXTERNAL DASUM, DLAMCH, DNRM2, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DBDSQR, DLASCL, XERBLA, ZGEBD2, ZLASCL, ZLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX, MIN

      ZQRT12 = ZERO;

      // Test that enough workspace is supplied

      if ( LWORK < M*N+2*min( M, N )+max( M, N ) ) {
         xerbla('ZQRT12', 7 );
         return;
      }

      // Quick return if possible

      MN = min( M, N );
      if (MN <= ZERO) return;

      NRMSVL = dnrm2( MN, S, 1 );

      // Copy upper triangle of A into work

      zlaset('Full', M, N, DCMPLX( ZERO ), DCMPLX( ZERO ), WORK, M );
      for (J = 1; J <= N; J++) {
         for (I = 1; I <= min( J, M ); I++) {
            WORK[( J-1 )*M+I] = A( I, J );
         }
      }

      // Get machine parameters

      SMLNUM = dlamch( 'S' ) / dlamch( 'P' );
      BIGNUM = ONE / SMLNUM;

      // Scale work if max entry outside range [SMLNUM,BIGNUM]

      ANRM = ZLANGE( 'M', M, N, WORK, M, DUMMY );
      ISCL = 0;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         zlascl('G', 0, 0, ANRM, SMLNUM, M, N, WORK, M, INFO );
         ISCL = 1;
      } else if ( ANRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         zlascl('G', 0, 0, ANRM, BIGNUM, M, N, WORK, M, INFO );
         ISCL = 1;
      }

      if ( ANRM != ZERO ) {

         // Compute SVD of work

         zgebd2(M, N, WORK, M, RWORK( 1 ), RWORK( MN+1 ), WORK( M*N+1 ), WORK( M*N+MN+1 ), WORK( M*N+2*MN+1 ), INFO );
         dbdsqr('Upper', MN, 0, 0, 0, RWORK( 1 ), RWORK( MN+1 ), DUMMY, MN, DUMMY, 1, DUMMY, MN, RWORK( 2*MN+1 ), INFO );

         if ( ISCL == 1 ) {
            if ( ANRM > BIGNUM ) {
               dlascl('G', 0, 0, BIGNUM, ANRM, MN, 1, RWORK( 1 ), MN, INFO );
            }
            if ( ANRM < SMLNUM ) {
               dlascl('G', 0, 0, SMLNUM, ANRM, MN, 1, RWORK( 1 ), MN, INFO );
            }
         }

      } else {

         for (I = 1; I <= MN; I++) {
            RWORK[I] = ZERO;
         }
      }

      // Compare s and singular values of work

      daxpy(MN, -ONE, S, 1, RWORK( 1 ), 1 );
      ZQRT12 = dasum( MN, RWORK( 1 ), 1 ) / ( dlamch( 'Epsilon' )*(max( M, N )).toDouble() );

      if (NRMSVL != ZERO) ZQRT12 = ZQRT12 / NRMSVL;

      }
