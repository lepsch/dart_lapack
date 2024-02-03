      SUBROUTINE CQRT15( SCALE, RKSEL, M, N, NRHS, A, LDA, B, LDB, S, RANK, NORMA, NORMB, ISEED, WORK, LWORK )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, LWORK, M, N, NRHS, RANK, RKSEL, SCALE;
      REAL               NORMA, NORMB
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      REAL               S( * )
      COMPLEX            A( LDA, * ), B( LDB, * ), WORK( LWORK )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO, SVMIN
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0, SVMIN = 0.1 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                INFO, J, MN;
      REAL               BIGNUM, EPS, SMLNUM, TEMP
      // ..
      // .. Local Arrays ..
      REAL               DUMMY( 1 )
      // ..
      // .. External Functions ..
      REAL               CLANGE, SASUM, SCNRM2, SLAMCH, SLARND
      // EXTERNAL CLANGE, SASUM, SCNRM2, SLAMCH, SLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CLARF, CLARNV, CLAROR, CLASCL, CLASET, CSSCAL, SLAORD, SLASCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, MAX, MIN
      // ..
      // .. Executable Statements ..

      MN = MIN( M, N )
      if ( LWORK < MAX( M+MN, MN*NRHS, 2*N+M ) ) {
         xerbla('CQRT15', 16 );
         RETURN
      }

      SMLNUM = SLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
      EPS = SLAMCH( 'Epsilon' )
      SMLNUM = ( SMLNUM / EPS ) / EPS
      BIGNUM = ONE / SMLNUM

      // Determine rank and (unscaled) singular values

      if ( RKSEL == 1 ) {
         RANK = MN
      } else if ( RKSEL == 2 ) {
         RANK = ( 3*MN ) / 4
         for (J = RANK + 1; J <= MN; J++) { // 10
            S( J ) = ZERO
         } // 10
      } else {
         xerbla('CQRT15', 2 );
      }

      if ( RANK > 0 ) {

         // Nontrivial case

         S( 1 ) = ONE
         for (J = 2; J <= RANK; J++) { // 30
            } // 20
            TEMP = SLARND( 1, ISEED )
            if ( TEMP > SVMIN ) {
               S( J ) = ABS( TEMP )
            } else {
               GO TO 20
            }
         } // 30
         slaord('Decreasing', RANK, S, 1 );

         // Generate 'rank' columns of a random orthogonal matrix in A

         clarnv(2, ISEED, M, WORK );
         csscal(M, ONE / SCNRM2( M, WORK, 1 ), WORK, 1 );
         claset('Full', M, RANK, CZERO, CONE, A, LDA );
         clarf('Left', M, RANK, WORK, 1, CMPLX( TWO ), A, LDA, WORK( M+1 ) );

         // workspace used: m+mn

         // Generate consistent rhs in the range space of A

         clarnv(2, ISEED, RANK*NRHS, WORK );
         cgemm('No transpose', 'No transpose', M, NRHS, RANK, CONE, A, LDA, WORK, RANK, CZERO, B, LDB );

         // work space used: <= mn *nrhs

         // generate (unscaled) matrix A

         for (J = 1; J <= RANK; J++) { // 40
            csscal(M, S( J ), A( 1, J ), 1 );
         } // 40
         if (RANK < N) CALL CLASET( 'Full', M, N-RANK, CZERO, CZERO, A( 1, RANK+1 ), LDA );
         claror('Right', 'No initialization', M, N, A, LDA, ISEED, WORK, INFO );

      } else {

         // work space used 2*n+m

         // Generate null matrix and rhs

         for (J = 1; J <= MN; J++) { // 50
            S( J ) = ZERO
         } // 50
         claset('Full', M, N, CZERO, CZERO, A, LDA );
         claset('Full', M, NRHS, CZERO, CZERO, B, LDB );

      }

      // Scale the matrix

      if ( SCALE != 1 ) {
         NORMA = CLANGE( 'Max', M, N, A, LDA, DUMMY )
         if ( NORMA != ZERO ) {
            if ( SCALE == 2 ) {

               // matrix scaled up

               clascl('General', 0, 0, NORMA, BIGNUM, M, N, A, LDA, INFO );
               slascl('General', 0, 0, NORMA, BIGNUM, MN, 1, S, MN, INFO );
               clascl('General', 0, 0, NORMA, BIGNUM, M, NRHS, B, LDB, INFO );
            } else if ( SCALE == 3 ) {

               // matrix scaled down

               clascl('General', 0, 0, NORMA, SMLNUM, M, N, A, LDA, INFO );
               slascl('General', 0, 0, NORMA, SMLNUM, MN, 1, S, MN, INFO );
               clascl('General', 0, 0, NORMA, SMLNUM, M, NRHS, B, LDB, INFO );
            } else {
               xerbla('CQRT15', 1 );
               RETURN
            }
         }
      }

      NORMA = SASUM( MN, S, 1 )
      NORMB = CLANGE( 'One-norm', M, NRHS, B, LDB, DUMMY )

      RETURN

      // End of CQRT15

      }
