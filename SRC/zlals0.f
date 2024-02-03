      SUBROUTINE ZLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B, LDB, BX, LDBX, PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, POLES, DIFL, DIFR, Z, K, C, S, RWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                GIVPTR, ICOMPQ, INFO, K, LDB, LDBX, LDGCOL, LDGNUM, NL, NR, NRHS, SQRE;
      double             C, S;
      // ..
      // .. Array Arguments ..
      int                GIVCOL( LDGCOL, * ), PERM( * );
      double             DIFL( * ), DIFR( LDGNUM, * ), GIVNUM( LDGNUM, * ), POLES( LDGNUM, * ), RWORK( * ), Z( * );
      COMPLEX*16         B( LDB, * ), BX( LDBX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO, NEGONE;
      const              ONE = 1.0D0, ZERO = 0.0D0, NEGONE = -1.0D0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, JCOL, JROW, M, N, NLP1;
      double             DIFLJ, DIFRJ, DJ, DSIGJ, DSIGJP, TEMP;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMV, XERBLA, ZCOPY, ZDROT, ZDSCAL, ZLACPY, ZLASCL
      // ..
      // .. External Functions ..
      double             DLAMC3, DNRM2;
      // EXTERNAL DLAMC3, DNRM2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, DIMAG, MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      N = NL + NR + 1

      if ( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) {
         INFO = -1
      } else if ( NL.LT.1 ) {
         INFO = -2
      } else if ( NR.LT.1 ) {
         INFO = -3
      } else if ( ( SQRE.LT.0 ) .OR. ( SQRE.GT.1 ) ) {
         INFO = -4
      } else if ( NRHS.LT.1 ) {
         INFO = -5
      } else if ( LDB.LT.N ) {
         INFO = -7
      } else if ( LDBX.LT.N ) {
         INFO = -9
      } else if ( GIVPTR.LT.0 ) {
         INFO = -11
      } else if ( LDGCOL.LT.N ) {
         INFO = -13
      } else if ( LDGNUM.LT.N ) {
         INFO = -15
      } else if ( K.LT.1 ) {
         INFO = -20
      }
      if ( INFO.NE.0 ) {
         xerbla('ZLALS0', -INFO );
         RETURN
      }

      M = N + SQRE
      NLP1 = NL + 1

      if ( ICOMPQ == 0 ) {

         // Apply back orthogonal transformations from the left.

         // Step (1L): apply back the Givens rotations performed.

         for (I = 1; I <= GIVPTR; I++) { // 10
            zdrot(NRHS, B( GIVCOL( I, 2 ), 1 ), LDB, B( GIVCOL( I, 1 ), 1 ), LDB, GIVNUM( I, 2 ), GIVNUM( I, 1 ) );
         } // 10

         // Step (2L): permute rows of B.

         zcopy(NRHS, B( NLP1, 1 ), LDB, BX( 1, 1 ), LDBX );
         for (I = 2; I <= N; I++) { // 20
            zcopy(NRHS, B( PERM( I ), 1 ), LDB, BX( I, 1 ), LDBX );
         } // 20

         // Step (3L): apply the inverse of the left singular vector
         // matrix to BX.

         if ( K == 1 ) {
            zcopy(NRHS, BX, LDBX, B, LDB );
            if ( Z( 1 ).LT.ZERO ) {
               zdscal(NRHS, NEGONE, B, LDB );
            }
         } else {
            for (J = 1; J <= K; J++) { // 100
               DIFLJ = DIFL( J )
               DJ = POLES( J, 1 )
               DSIGJ = -POLES( J, 2 )
               if ( J.LT.K ) {
                  DIFRJ = -DIFR( J, 1 )
                  DSIGJP = -POLES( J+1, 2 )
               }
               if ( ( Z( J ) == ZERO ) .OR. ( POLES( J, 2 ) == ZERO ) ) {
                  RWORK( J ) = ZERO
               } else {
                  RWORK( J ) = -POLES( J, 2 )*Z( J ) / DIFLJ / ( POLES( J, 2 )+DJ )
               }
               for (I = 1; I <= J - 1; I++) { // 30
                  if ( ( Z( I ) == ZERO ) .OR. ( POLES( I, 2 ) == ZERO ) ) {
                     RWORK( I ) = ZERO
                  } else {

                     // Use calls to the subroutine DLAMC3 to enforce the
                     // parentheses (x+y)+z. The goal is to prevent
                     // optimizing compilers from doing x+(y+z).

                     RWORK( I ) = POLES( I, 2 )*Z( I ) / ( DLAMC3( POLES( I, 2 ), DSIGJ )- DIFLJ ) / ( POLES( I, 2 )+DJ )
                  }
               } // 30
               for (I = J + 1; I <= K; I++) { // 40
                  if ( ( Z( I ) == ZERO ) .OR. ( POLES( I, 2 ) == ZERO ) ) {
                     RWORK( I ) = ZERO
                  } else {
                     RWORK( I ) = POLES( I, 2 )*Z( I ) / ( DLAMC3( POLES( I, 2 ), DSIGJP )+ DIFRJ ) / ( POLES( I, 2 )+DJ )
                  }
               } // 40
               RWORK( 1 ) = NEGONE
               TEMP = DNRM2( K, RWORK, 1 )

               // Since B and BX are complex, the following call to DGEMV
               // is performed in two steps (real and imaginary parts).

               // CALL DGEMV( 'T', K, NRHS, ONE, BX, LDBX, WORK, 1, ZERO,
*    $                     B( J, 1 ), LDB )

               I = K + NRHS*2
               for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 60
                  for (JROW = 1; JROW <= K; JROW++) { // 50
                     I = I + 1
                     RWORK( I ) = DBLE( BX( JROW, JCOL ) )
                  } // 50
               } // 60
               dgemv('T', K, NRHS, ONE, RWORK( 1+K+NRHS*2 ), K, RWORK( 1 ), 1, ZERO, RWORK( 1+K ), 1 );
               I = K + NRHS*2
               for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 80
                  for (JROW = 1; JROW <= K; JROW++) { // 70
                     I = I + 1
                     RWORK( I ) = DIMAG( BX( JROW, JCOL ) )
                  } // 70
               } // 80
               dgemv('T', K, NRHS, ONE, RWORK( 1+K+NRHS*2 ), K, RWORK( 1 ), 1, ZERO, RWORK( 1+K+NRHS ), 1 );
               for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 90
                  B( J, JCOL ) = DCMPLX( RWORK( JCOL+K ), RWORK( JCOL+K+NRHS ) )
               } // 90
               zlascl('G', 0, 0, TEMP, ONE, 1, NRHS, B( J, 1 ), LDB, INFO );
            } // 100
         }

         // Move the deflated rows of BX to B also.

         IF( K.LT.MAX( M, N ) ) CALL ZLACPY( 'A', N-K, NRHS, BX( K+1, 1 ), LDBX, B( K+1, 1 ), LDB )
      } else {

         // Apply back the right orthogonal transformations.

         // Step (1R): apply back the new right singular vector matrix
         // to B.

         if ( K == 1 ) {
            zcopy(NRHS, B, LDB, BX, LDBX );
         } else {
            for (J = 1; J <= K; J++) { // 180
               DSIGJ = POLES( J, 2 )
               if ( Z( J ) == ZERO ) {
                  RWORK( J ) = ZERO
               } else {
                  RWORK( J ) = -Z( J ) / DIFL( J ) / ( DSIGJ+POLES( J, 1 ) ) / DIFR( J, 2 )
               }
               for (I = 1; I <= J - 1; I++) { // 110
                  if ( Z( J ) == ZERO ) {
                     RWORK( I ) = ZERO
                  } else {

                     // Use calls to the subroutine DLAMC3 to enforce the
                     // parentheses (x+y)+z. The goal is to prevent
                     // optimizing compilers from doing x+(y+z).

                     RWORK( I ) = Z( J ) / ( DLAMC3( DSIGJ, -POLES( I+1, 2 ) )-DIFR( I, 1 ) ) / ( DSIGJ+POLES( I, 1 ) ) / DIFR( I, 2 )
                  }
               } // 110
               for (I = J + 1; I <= K; I++) { // 120
                  if ( Z( J ) == ZERO ) {
                     RWORK( I ) = ZERO
                  } else {
                     RWORK( I ) = Z( J ) / ( DLAMC3( DSIGJ, -POLES( I, 2 ) )-DIFL( I ) ) / ( DSIGJ+POLES( I, 1 ) ) / DIFR( I, 2 )
                  }
               } // 120

               // Since B and BX are complex, the following call to DGEMV
               // is performed in two steps (real and imaginary parts).

               // CALL DGEMV( 'T', K, NRHS, ONE, B, LDB, WORK, 1, ZERO,
*    $                     BX( J, 1 ), LDBX )

               I = K + NRHS*2
               for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 140
                  for (JROW = 1; JROW <= K; JROW++) { // 130
                     I = I + 1
                     RWORK( I ) = DBLE( B( JROW, JCOL ) )
                  } // 130
               } // 140
               dgemv('T', K, NRHS, ONE, RWORK( 1+K+NRHS*2 ), K, RWORK( 1 ), 1, ZERO, RWORK( 1+K ), 1 );
               I = K + NRHS*2
               for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 160
                  for (JROW = 1; JROW <= K; JROW++) { // 150
                     I = I + 1
                     RWORK( I ) = DIMAG( B( JROW, JCOL ) )
                  } // 150
               } // 160
               dgemv('T', K, NRHS, ONE, RWORK( 1+K+NRHS*2 ), K, RWORK( 1 ), 1, ZERO, RWORK( 1+K+NRHS ), 1 );
               for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 170
                  BX( J, JCOL ) = DCMPLX( RWORK( JCOL+K ), RWORK( JCOL+K+NRHS ) )
               } // 170
            } // 180
         }

         // Step (2R): if SQRE = 1, apply back the rotation that is
         // related to the right null space of the subproblem.

         if ( SQRE == 1 ) {
            zcopy(NRHS, B( M, 1 ), LDB, BX( M, 1 ), LDBX );
            zdrot(NRHS, BX( 1, 1 ), LDBX, BX( M, 1 ), LDBX, C, S );
         }
         IF( K.LT.MAX( M, N ) ) CALL ZLACPY( 'A', N-K, NRHS, B( K+1, 1 ), LDB, BX( K+1, 1 ), LDBX )

         // Step (3R): permute rows of B.

         zcopy(NRHS, BX( 1, 1 ), LDBX, B( NLP1, 1 ), LDB );
         if ( SQRE == 1 ) {
            zcopy(NRHS, BX( M, 1 ), LDBX, B( M, 1 ), LDB );
         }
         for (I = 2; I <= N; I++) { // 190
            zcopy(NRHS, BX( I, 1 ), LDBX, B( PERM( I ), 1 ), LDB );
         } // 190

         // Step (4R): apply back the Givens rotations performed.

         DO 200 I = GIVPTR, 1, -1
            zdrot(NRHS, B( GIVCOL( I, 2 ), 1 ), LDB, B( GIVCOL( I, 1 ), 1 ), LDB, GIVNUM( I, 2 ), -GIVNUM( I, 1 ) );
         } // 200
      }

      RETURN

      // End of ZLALS0

      }
