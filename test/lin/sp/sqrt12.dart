      double sqrt12(final int M, final int N, final Matrix<double> A, final int LDA, final int S, final Array<double> WORK, final int LWORK,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LWORK, M, N;
      double               A( LDA, * ), S( * ), WORK( LWORK );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, INFO, ISCL, J, MN;
      double               ANRM, BIGNUM, NRMSVL, SMLNUM;
      // ..
      // .. External Functions ..
      //- REAL               SASUM, SLAMCH, SLANGE, SNRM2;
      // EXTERNAL SASUM, SLAMCH, SLANGE, SNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SBDSQR, SGEBD2, SLASCL, SLASET, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      double               DUMMY( 1 );

      SQRT12 = ZERO;

      // Test that enough workspace is supplied

      if ( LWORK < max( M*N+4*min( M, N )+max( M, N ), M*N+2*min( M, N )+4*N) ) {
         xerbla('SQRT12', 7 );
         return;
      }

      // Quick return if possible

      MN = min( M, N );
      if (MN <= ZERO) return;

      NRMSVL = SNRM2( MN, S, 1 );

      // Copy upper triangle of A into work

      slaset('Full', M, N, ZERO, ZERO, WORK, M );
      for (J = 1; J <= N; J++) {
         for (I = 1; I <= min( J, M ); I++) {
            WORK[( J-1 )*M+I] = A( I, J );
         }
      }

      // Get machine parameters

      SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'P' );
      BIGNUM = ONE / SMLNUM;

      // Scale work if max entry outside range [SMLNUM,BIGNUM]

      ANRM = SLANGE( 'M', M, N, WORK, M, DUMMY );
      ISCL = 0;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         slascl('G', 0, 0, ANRM, SMLNUM, M, N, WORK, M, INFO );
         ISCL = 1;
      } else if ( ANRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         slascl('G', 0, 0, ANRM, BIGNUM, M, N, WORK, M, INFO );
         ISCL = 1;
      }

      if ( ANRM != ZERO ) {

         // Compute SVD of work

         sgebd2(M, N, WORK, M, WORK( M*N+1 ), WORK( M*N+MN+1 ), WORK( M*N+2*MN+1 ), WORK( M*N+3*MN+1 ), WORK( M*N+4*MN+1 ), INFO );
         sbdsqr('Upper', MN, 0, 0, 0, WORK( M*N+1 ), WORK( M*N+MN+1 ), DUMMY, MN, DUMMY, 1, DUMMY, MN, WORK( M*N+2*MN+1 ), INFO );

         if ( ISCL == 1 ) {
            if ( ANRM > BIGNUM ) {
               slascl('G', 0, 0, BIGNUM, ANRM, MN, 1, WORK( M*N+1 ), MN, INFO );
            }
            if ( ANRM < SMLNUM ) {
               slascl('G', 0, 0, SMLNUM, ANRM, MN, 1, WORK( M*N+1 ), MN, INFO );
            }
         }

      } else {

         for (I = 1; I <= MN; I++) {
            WORK[M*N+I] = ZERO;
         }
      }

      // Compare s and singular values of work

      saxpy(MN, -ONE, S, 1, WORK( M*N+1 ), 1 );
      SQRT12 = SASUM( MN, WORK( M*N+1 ), 1 ) / ( SLAMCH( 'Epsilon' )*REAL( max( M, N ) ) )       IF( NRMSVL != ZERO ) SQRT12 = SQRT12 / NRMSVL;

      }
