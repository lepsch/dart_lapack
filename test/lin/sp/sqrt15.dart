      void sqrt15(SCALE, RKSEL, M, N, NRHS, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, S, RANK, NORMA, NORMB, final Array<int> ISEED, final Array<double> WORK, final int LWORK) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDB, LWORK, M, N, NRHS, RANK, RKSEL, SCALE;
      double               NORMA, NORMB;
      int                ISEED( 4 );
      double               A( LDA, * ), B( LDB, * ), S( * ), WORK( LWORK );
      // ..

      double               ZERO, ONE, TWO, SVMIN;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0, SVMIN = 0.1 ;
      int                INFO, J, MN;
      double               BIGNUM, EPS, SMLNUM, TEMP;
      double               DUMMY( 1 );
      // ..
      // .. External Functions ..
      //- REAL               SASUM, SLAMCH, SLANGE, SLARND, SNRM2;
      // EXTERNAL SASUM, SLAMCH, SLANGE, SLARND, SNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLAORD, SLARF, SLARNV, SLAROR, SLASCL, SLASET, SSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN

      MN = min( M, N );
      if ( LWORK < max( M+MN, MN*NRHS, 2*N+M ) ) {
         xerbla('SQRT15', 16 );
         return;
      }

      SMLNUM = SLAMCH( 'Safe minimum' );
      BIGNUM = ONE / SMLNUM;
      EPS = SLAMCH( 'Epsilon' );
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
         xerbla('SQRT15', 2 );
      }

      if ( RANK > 0 ) {

         // Nontrivial case

         S[1] = ONE;
         for (J = 2; J <= RANK; J++) { // 30
            } // 20
            TEMP = SLARND( 1, ISEED );
            if ( TEMP > SVMIN ) {
               S[J] = ( TEMP ).abs();
            } else {
               GO TO 20;
            }
         } // 30
         slaord('Decreasing', RANK, S, 1 );

         // Generate 'rank' columns of a random orthogonal matrix in A

         slarnv(2, ISEED, M, WORK );
         sscal(M, ONE / SNRM2( M, WORK, 1 ), WORK, 1 );
         slaset('Full', M, RANK, ZERO, ONE, A, LDA );
         slarf('Left', M, RANK, WORK, 1, TWO, A, LDA, WORK( M+1 ) );

         // workspace used: m+mn

         // Generate consistent rhs in the range space of A

         slarnv(2, ISEED, RANK*NRHS, WORK );
         sgemm('No transpose', 'No transpose', M, NRHS, RANK, ONE, A, LDA, WORK, RANK, ZERO, B, LDB );

         // work space used: <= mn *nrhs

         // generate (unscaled) matrix A

         for (J = 1; J <= RANK; J++) { // 40
            sscal(M, S( J ), A( 1, J ), 1 );
         } // 40
         if (RANK < N) slaset( 'Full', M, N-RANK, ZERO, ZERO, A( 1, RANK+1 ), LDA );
         slaror('Right', 'No initialization', M, N, A, LDA, ISEED, WORK, INFO );

      } else {

         // work space used 2*n+m

         // Generate null matrix and rhs

         for (J = 1; J <= MN; J++) { // 50
            S[J] = ZERO;
         } // 50
         slaset('Full', M, N, ZERO, ZERO, A, LDA );
         slaset('Full', M, NRHS, ZERO, ZERO, B, LDB );

      }

      // Scale the matrix

      if ( SCALE != 1 ) {
         NORMA = SLANGE( 'Max', M, N, A, LDA, DUMMY );
         if ( NORMA != ZERO ) {
            if ( SCALE == 2 ) {

               // matrix scaled up

               slascl('General', 0, 0, NORMA, BIGNUM, M, N, A, LDA, INFO );
               slascl('General', 0, 0, NORMA, BIGNUM, MN, 1, S, MN, INFO );
               slascl('General', 0, 0, NORMA, BIGNUM, M, NRHS, B, LDB, INFO );
            } else if ( SCALE == 3 ) {

               // matrix scaled down

               slascl('General', 0, 0, NORMA, SMLNUM, M, N, A, LDA, INFO );
               slascl('General', 0, 0, NORMA, SMLNUM, MN, 1, S, MN, INFO );
               slascl('General', 0, 0, NORMA, SMLNUM, M, NRHS, B, LDB, INFO );
            } else {
               xerbla('SQRT15', 1 );
               return;
            }
         }
      }

      NORMA = SASUM( MN, S, 1 );
      NORMB = SLANGE( 'One-norm', M, NRHS, B, LDB, DUMMY );

      }
