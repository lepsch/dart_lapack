      void zqrt15(final int SCALE, final int RKSEL, final int M, final int N, final int NRHS, final Matrix<double> A_, final int LDA, final Matrix<double> B_, final int LDB, final int S, final int RANK, final int NORMA, final int NORMB, final Array<int> ISEED_, final Array<double> WORK_, final int LWORK,) {
  final A = A_.dim();
  final B = B_.dim();
  final ISEED = ISEED_.dim();
  final WORK = WORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDB, LWORK, M, N, NRHS, RANK, RKSEL, SCALE;
      double             NORMA, NORMB;
      int                ISEED( 4 );
      double             S( * );
      Complex         A( LDA, * ), B( LDB, * ), WORK( LWORK );
      // ..

      double             ZERO, ONE, TWO, SVMIN;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0, SVMIN = 0.1 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                INFO, J, MN;
      double             BIGNUM, EPS, SMLNUM, TEMP;
      double             DUMMY( 1 );
      // ..
      // .. External Functions ..
      //- double             DASUM, DLAMCH, DLARND, DZNRM2, ZLANGE;
      // EXTERNAL DASUM, DLAMCH, DLARND, DZNRM2, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLAORD, DLASCL, XERBLA, ZDSCAL, ZGEMM, ZLARF, ZLARNV, ZLAROR, ZLASCL, ZLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCMPLX, MAX, MIN

      MN = min( M, N );
      if ( LWORK < max( M+MN, MN*NRHS, 2*N+M ) ) {
         xerbla('ZQRT15', 16 );
         return;
      }

      SMLNUM = dlamch( 'Safe minimum' );
      BIGNUM = ONE / SMLNUM;
      EPS = dlamch( 'Epsilon' );
      SMLNUM = ( SMLNUM / EPS ) / EPS;
      BIGNUM = ONE / SMLNUM;

      // Determine rank and (unscaled) singular values

      if ( RKSEL == 1 ) {
         RANK = MN;
      } else if ( RKSEL == 2 ) {
         RANK = ( 3*MN ) / 4;
         for (J = RANK + 1; J <= MN; J++) { // 10
            S[J] = ZERO;
         } // 10
      } else {
         xerbla('ZQRT15', 2 );
      }

      if ( RANK > 0 ) {

         // Nontrivial case

         S[1] = ONE;
         for (J = 2; J <= RANK; J++) { // 30
            } // 20
            TEMP = dlarnd( 1, ISEED );
            if ( TEMP > SVMIN ) {
               S[J] = ( TEMP ).abs();
            } else {
               GO TO 20;
            }
         } // 30
         dlaord('Decreasing', RANK, S, 1 );

         // Generate 'rank' columns of a random orthogonal matrix in A

         zlarnv(2, ISEED, M, WORK );
         zdscal(M, ONE / DZNRM2( M, WORK, 1 ), WORK, 1 );
         zlaset('Full', M, RANK, CZERO, CONE, A, LDA );
         zlarf('Left', M, RANK, WORK, 1, DCMPLX( TWO ), A, LDA, WORK( M+1 ) );

         // workspace used: m+mn

         // Generate consistent rhs in the range space of A

         zlarnv(2, ISEED, RANK*NRHS, WORK );
         zgemm('No transpose', 'No transpose', M, NRHS, RANK, CONE, A, LDA, WORK, RANK, CZERO, B, LDB );

         // work space used: <= mn *nrhs

         // generate (unscaled) matrix A

         for (J = 1; J <= RANK; J++) { // 40
            zdscal(M, S( J ), A( 1, J ), 1 );
         } // 40
         if (RANK < N) zlaset( 'Full', M, N-RANK, CZERO, CZERO, A( 1, RANK+1 ), LDA );
         zlaror('Right', 'No initialization', M, N, A, LDA, ISEED, WORK, INFO );

      } else {

         // work space used 2*n+m

         // Generate null matrix and rhs

         for (J = 1; J <= MN; J++) { // 50
            S[J] = ZERO;
         } // 50
         zlaset('Full', M, N, CZERO, CZERO, A, LDA );
         zlaset('Full', M, NRHS, CZERO, CZERO, B, LDB );

      }

      // Scale the matrix

      if ( SCALE != 1 ) {
         NORMA = ZLANGE( 'Max', M, N, A, LDA, DUMMY );
         if ( NORMA != ZERO ) {
            if ( SCALE == 2 ) {

               // matrix scaled up

               zlascl('General', 0, 0, NORMA, BIGNUM, M, N, A, LDA, INFO );
               dlascl('General', 0, 0, NORMA, BIGNUM, MN, 1, S, MN, INFO );
               zlascl('General', 0, 0, NORMA, BIGNUM, M, NRHS, B, LDB, INFO );
            } else if ( SCALE == 3 ) {

               // matrix scaled down

               zlascl('General', 0, 0, NORMA, SMLNUM, M, N, A, LDA, INFO );
               dlascl('General', 0, 0, NORMA, SMLNUM, MN, 1, S, MN, INFO );
               zlascl('General', 0, 0, NORMA, SMLNUM, M, NRHS, B, LDB, INFO );
            } else {
               xerbla('ZQRT15', 1 );
               return;
            }
         }
      }

      NORMA = dasum( MN, S, 1 );
      NORMB = ZLANGE( 'One-norm', M, NRHS, B, LDB, DUMMY );

      }
