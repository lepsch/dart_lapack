      SUBROUTINE DQRT15( SCALE, RKSEL, M, N, NRHS, A, LDA, B, LDB, S, RANK, NORMA, NORMB, ISEED, WORK, LWORK )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, LWORK, M, N, NRHS, RANK, RKSEL, SCALE;
      double             NORMA, NORMB;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      double             A( LDA, * ), B( LDB, * ), S( * ), WORK( LWORK );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO, SVMIN;
      const              ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, SVMIN = 0.1D0 ;
      // ..
      // .. Local Scalars ..
      int                INFO, J, MN;
      double             BIGNUM, EPS, SMLNUM, TEMP;
      // ..
      // .. Local Arrays ..
      double             DUMMY( 1 );
      // ..
      // .. External Functions ..
      double             DASUM, DLAMCH, DLANGE, DLARND, DNRM2;
      // EXTERNAL DASUM, DLAMCH, DLANGE, DLARND, DNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLAORD, DLARF, DLARNV, DLAROR, DLASCL, DLASET, DSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      MN = MIN( M, N )
      if ( LWORK.LT.MAX( M+MN, MN*NRHS, 2*N+M ) ) {
         xerbla('DQRT15', 16 );
         RETURN
      }

      SMLNUM = DLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
      EPS = DLAMCH( 'Epsilon' )
      SMLNUM = ( SMLNUM / EPS ) / EPS
      BIGNUM = ONE / SMLNUM

      // Determine rank and (unscaled) singular values

      if ( RKSEL.EQ.1 ) {
         RANK = MN
      } else if ( RKSEL.EQ.2 ) {
         RANK = ( 3*MN ) / 4
         DO 10 J = RANK + 1, MN
            S( J ) = ZERO
   10    CONTINUE
      } else {
         xerbla('DQRT15', 2 );
      }

      if ( RANK.GT.0 ) {

         // Nontrivial case

         S( 1 ) = ONE
         for (J = 2; J <= RANK; J++) { // 30
   20       CONTINUE
            TEMP = DLARND( 1, ISEED )
            if ( TEMP.GT.SVMIN ) {
               S( J ) = ABS( TEMP )
            } else {
               GO TO 20
            }
   30    CONTINUE
         dlaord('Decreasing', RANK, S, 1 );

         // Generate 'rank' columns of a random orthogonal matrix in A

         dlarnv(2, ISEED, M, WORK );
         dscal(M, ONE / DNRM2( M, WORK, 1 ), WORK, 1 );
         dlaset('Full', M, RANK, ZERO, ONE, A, LDA );
         dlarf('Left', M, RANK, WORK, 1, TWO, A, LDA, WORK( M+1 ) );

         // workspace used: m+mn

         // Generate consistent rhs in the range space of A

         dlarnv(2, ISEED, RANK*NRHS, WORK );
         dgemm('No transpose', 'No transpose', M, NRHS, RANK, ONE, A, LDA, WORK, RANK, ZERO, B, LDB );

         // work space used: <= mn *nrhs

         // generate (unscaled) matrix A

         for (J = 1; J <= RANK; J++) { // 40
            dscal(M, S( J ), A( 1, J ), 1 );
   40    CONTINUE
         IF( RANK.LT.N ) CALL DLASET( 'Full', M, N-RANK, ZERO, ZERO, A( 1, RANK+1 ), LDA )
         dlaror('Right', 'No initialization', M, N, A, LDA, ISEED, WORK, INFO );

      } else {

         // work space used 2*n+m

         // Generate null matrix and rhs

         for (J = 1; J <= MN; J++) { // 50
            S( J ) = ZERO
   50    CONTINUE
         dlaset('Full', M, N, ZERO, ZERO, A, LDA );
         dlaset('Full', M, NRHS, ZERO, ZERO, B, LDB );

      }

      // Scale the matrix

      if ( SCALE.NE.1 ) {
         NORMA = DLANGE( 'Max', M, N, A, LDA, DUMMY )
         if ( NORMA.NE.ZERO ) {
            if ( SCALE.EQ.2 ) {

               // matrix scaled up

               dlascl('General', 0, 0, NORMA, BIGNUM, M, N, A, LDA, INFO )                CALL DLASCL( 'General', 0, 0, NORMA, BIGNUM, MN, 1, S, MN, INFO )                CALL DLASCL( 'General', 0, 0, NORMA, BIGNUM, M, NRHS, B, LDB, INFO );
            } else if ( SCALE.EQ.3 ) {

               // matrix scaled down

               dlascl('General', 0, 0, NORMA, SMLNUM, M, N, A, LDA, INFO )                CALL DLASCL( 'General', 0, 0, NORMA, SMLNUM, MN, 1, S, MN, INFO )                CALL DLASCL( 'General', 0, 0, NORMA, SMLNUM, M, NRHS, B, LDB, INFO );
            } else {
               xerbla('DQRT15', 1 );
               RETURN
            }
         }
      }

      NORMA = DASUM( MN, S, 1 )
      NORMB = DLANGE( 'One-norm', M, NRHS, B, LDB, DUMMY )

      RETURN

      // End of DQRT15

      }
