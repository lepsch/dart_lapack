      REAL             FUNCTION CQRT12( M, N, A, LDA, S, WORK, LWORK, RWORK );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * ), S( * );
      COMPLEX            A( LDA, * ), WORK( LWORK );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, ISCL, J, MN;
      REAL               ANRM, BIGNUM, NRMSVL, SMLNUM;
      // ..
      // .. Local Arrays ..
      REAL               DUMMY( 1 );
      // ..
      // .. External Functions ..
      REAL               CLANGE, SASUM, SLAMCH, SNRM2;
      // EXTERNAL CLANGE, SASUM, SLAMCH, SNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEBD2, CLASCL, CLASET, SAXPY, SBDSQR, SLASCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      CQRT12 = ZERO;

      // Test that enough workspace is supplied

      if ( LWORK < M*N+2*MIN( M, N )+MAX( M, N ) ) {
         xerbla('CQRT12', 7 );
         return;
      }

      // Quick return if possible

      MN = MIN( M, N );
      if (MN <= ZERO) RETURN;

      NRMSVL = SNRM2( MN, S, 1 );

      // Copy upper triangle of A into work

      claset('Full', M, N, CMPLX( ZERO ), CMPLX( ZERO ), WORK, M );
      for (J = 1; J <= N; J++) {
         DO I = 1, MIN( J, M );
            WORK( ( J-1 )*M+I ) = A( I, J );
         }
      }

      // Get machine parameters

      SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'P' );
      BIGNUM = ONE / SMLNUM;

      // Scale work if max entry outside range [SMLNUM,BIGNUM]

      ANRM = CLANGE( 'M', M, N, WORK, M, DUMMY );
      ISCL = 0;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         clascl('G', 0, 0, ANRM, SMLNUM, M, N, WORK, M, INFO );
         ISCL = 1;
      } else if ( ANRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         clascl('G', 0, 0, ANRM, BIGNUM, M, N, WORK, M, INFO );
         ISCL = 1;
      }

      if ( ANRM != ZERO ) {

         // Compute SVD of work

         cgebd2(M, N, WORK, M, RWORK( 1 ), RWORK( MN+1 ), WORK( M*N+1 ), WORK( M*N+MN+1 ), WORK( M*N+2*MN+1 ), INFO );
         sbdsqr('Upper', MN, 0, 0, 0, RWORK( 1 ), RWORK( MN+1 ), DUMMY, MN, DUMMY, 1, DUMMY, MN, RWORK( 2*MN+1 ), INFO );

         if ( ISCL == 1 ) {
            if ( ANRM > BIGNUM ) {
               slascl('G', 0, 0, BIGNUM, ANRM, MN, 1, RWORK( 1 ), MN, INFO );
            }
            if ( ANRM < SMLNUM ) {
               slascl('G', 0, 0, SMLNUM, ANRM, MN, 1, RWORK( 1 ), MN, INFO );
            }
         }

      } else {

         for (I = 1; I <= MN; I++) {
            RWORK( I ) = ZERO;
         }
      }

      // Compare s and singular values of work

      saxpy(MN, -ONE, S, 1, RWORK( 1 ), 1 );
      CQRT12 = SASUM( MN, RWORK( 1 ), 1 ) / ( SLAMCH( 'Epsilon' )*REAL( MAX( M, N ) ) )       IF( NRMSVL != ZERO ) CQRT12 = CQRT12 / NRMSVL;

      return;

      // End of CQRT12

      }
