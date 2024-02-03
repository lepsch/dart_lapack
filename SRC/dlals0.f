      void dlals0(ICOMPQ, NL, NR, SQRE, NRHS, B, LDB, BX, LDBX, PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, POLES, DIFL, DIFR, Z, K, C, S, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                GIVPTR, ICOMPQ, INFO, K, LDB, LDBX, LDGCOL, LDGNUM, NL, NR, NRHS, SQRE;
      double             C, S;
      // ..
      // .. Array Arguments ..
      int                GIVCOL( LDGCOL, * ), PERM( * );
      double             B( LDB, * ), BX( LDBX, * ), DIFL( * ), DIFR( LDGNUM, * ), GIVNUM( LDGNUM, * ), POLES( LDGNUM, * ), WORK( * ), Z( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO, NEGONE;
      const              ONE = 1.0, ZERO = 0.0, NEGONE = -1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, M, N, NLP1;
      double             DIFLJ, DIFRJ, DJ, DSIGJ, DSIGJP, TEMP;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGEMV, DLACPY, DLASCL, DROT, DSCAL, XERBLA
      // ..
      // .. External Functions ..
      double             DLAMC3, DNRM2;
      // EXTERNAL DLAMC3, DNRM2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      N = NL + NR + 1;

      if ( ( ICOMPQ < 0 ) || ( ICOMPQ > 1 ) ) {
         INFO = -1;
      } else if ( NL < 1 ) {
         INFO = -2;
      } else if ( NR < 1 ) {
         INFO = -3;
      } else if ( ( SQRE < 0 ) || ( SQRE > 1 ) ) {
         INFO = -4;
      } else if ( NRHS < 1 ) {
         INFO = -5;
      } else if ( LDB < N ) {
         INFO = -7;
      } else if ( LDBX < N ) {
         INFO = -9;
      } else if ( GIVPTR < 0 ) {
         INFO = -11;
      } else if ( LDGCOL < N ) {
         INFO = -13;
      } else if ( LDGNUM < N ) {
         INFO = -15;
      } else if ( K < 1 ) {
         INFO = -20;
      }
      if ( INFO != 0 ) {
         xerbla('DLALS0', -INFO );
         return;
      }

      M = N + SQRE;
      NLP1 = NL + 1;

      if ( ICOMPQ == 0 ) {

         // Apply back orthogonal transformations from the left.

         // Step (1L): apply back the Givens rotations performed.

         for (I = 1; I <= GIVPTR; I++) { // 10
            drot(NRHS, B( GIVCOL( I, 2 ), 1 ), LDB, B( GIVCOL( I, 1 ), 1 ), LDB, GIVNUM( I, 2 ), GIVNUM( I, 1 ) );
         } // 10

         // Step (2L): permute rows of B.

         dcopy(NRHS, B( NLP1, 1 ), LDB, BX( 1, 1 ), LDBX );
         for (I = 2; I <= N; I++) { // 20
            dcopy(NRHS, B( PERM( I ), 1 ), LDB, BX( I, 1 ), LDBX );
         } // 20

         // Step (3L): apply the inverse of the left singular vector
         // matrix to BX.

         if ( K == 1 ) {
            dcopy(NRHS, BX, LDBX, B, LDB );
            if ( Z( 1 ) < ZERO ) {
               dscal(NRHS, NEGONE, B, LDB );
            }
         } else {
            for (J = 1; J <= K; J++) { // 50
               DIFLJ = DIFL( J );
               DJ = POLES( J, 1 );
               DSIGJ = -POLES( J, 2 );
               if ( J < K ) {
                  DIFRJ = -DIFR( J, 1 );
                  DSIGJP = -POLES( J+1, 2 );
               }
               if ( ( Z( J ) == ZERO ) || ( POLES( J, 2 ) == ZERO ) ) {
                  WORK( J ) = ZERO;
               } else {
                  WORK( J ) = -POLES( J, 2 )*Z( J ) / DIFLJ / ( POLES( J, 2 )+DJ );
               }
               for (I = 1; I <= J - 1; I++) { // 30
                  if ( ( Z( I ) == ZERO ) || ( POLES( I, 2 ) == ZERO ) ) {
                     WORK( I ) = ZERO;
                  } else {

                     // Use calls to the subroutine DLAMC3 to enforce the
                     // parentheses (x+y)+z. The goal is to prevent
                     // optimizing compilers from doing x+(y+z).

                     WORK( I ) = POLES( I, 2 )*Z( I ) / ( DLAMC3( POLES( I, 2 ), DSIGJ )- DIFLJ ) / ( POLES( I, 2 )+DJ );
                  }
               } // 30
               for (I = J + 1; I <= K; I++) { // 40
                  if ( ( Z( I ) == ZERO ) || ( POLES( I, 2 ) == ZERO ) ) {
                     WORK( I ) = ZERO;
                  } else {
                     WORK( I ) = POLES( I, 2 )*Z( I ) / ( DLAMC3( POLES( I, 2 ), DSIGJP )+ DIFRJ ) / ( POLES( I, 2 )+DJ );
                  }
               } // 40
               WORK( 1 ) = NEGONE;
               TEMP = DNRM2( K, WORK, 1 );
               dgemv('T', K, NRHS, ONE, BX, LDBX, WORK, 1, ZERO, B( J, 1 ), LDB );
               dlascl('G', 0, 0, TEMP, ONE, 1, NRHS, B( J, 1 ), LDB, INFO );
            } // 50
         }

         // Move the deflated rows of BX to B also.

         if( K < max( M, N ) ) CALL DLACPY( 'A', N-K, NRHS, BX( K+1, 1 ), LDBX, B( K+1, 1 ), LDB );
      } else {

         // Apply back the right orthogonal transformations.

         // Step (1R): apply back the new right singular vector matrix
         // to B.

         if ( K == 1 ) {
            dcopy(NRHS, B, LDB, BX, LDBX );
         } else {
            for (J = 1; J <= K; J++) { // 80
               DSIGJ = POLES( J, 2 );
               if ( Z( J ) == ZERO ) {
                  WORK( J ) = ZERO;
               } else {
                  WORK( J ) = -Z( J ) / DIFL( J ) / ( DSIGJ+POLES( J, 1 ) ) / DIFR( J, 2 );
               }
               for (I = 1; I <= J - 1; I++) { // 60
                  if ( Z( J ) == ZERO ) {
                     WORK( I ) = ZERO;
                  } else {

                     // Use calls to the subroutine DLAMC3 to enforce the
                     // parentheses (x+y)+z. The goal is to prevent
                     // optimizing compilers from doing x+(y+z).

                     WORK( I ) = Z( J ) / ( DLAMC3( DSIGJ, -POLES( I+1, 2 ) )-DIFR( I, 1 ) ) / ( DSIGJ+POLES( I, 1 ) ) / DIFR( I, 2 );
                  }
               } // 60
               for (I = J + 1; I <= K; I++) { // 70
                  if ( Z( J ) == ZERO ) {
                     WORK( I ) = ZERO;
                  } else {
                     WORK( I ) = Z( J ) / ( DLAMC3( DSIGJ, -POLES( I, 2 ) )-DIFL( I ) ) / ( DSIGJ+POLES( I, 1 ) ) / DIFR( I, 2 );
                  }
               } // 70
               dgemv('T', K, NRHS, ONE, B, LDB, WORK, 1, ZERO, BX( J, 1 ), LDBX );
            } // 80
         }

         // Step (2R): if SQRE = 1, apply back the rotation that is
         // related to the right null space of the subproblem.

         if ( SQRE == 1 ) {
            dcopy(NRHS, B( M, 1 ), LDB, BX( M, 1 ), LDBX );
            drot(NRHS, BX( 1, 1 ), LDBX, BX( M, 1 ), LDBX, C, S );
         }
         if( K < max( M, N ) ) CALL DLACPY( 'A', N-K, NRHS, B( K+1, 1 ), LDB, BX( K+1, 1 ), LDBX );

         // Step (3R): permute rows of B.

         dcopy(NRHS, BX( 1, 1 ), LDBX, B( NLP1, 1 ), LDB );
         if ( SQRE == 1 ) {
            dcopy(NRHS, BX( M, 1 ), LDBX, B( M, 1 ), LDB );
         }
         for (I = 2; I <= N; I++) { // 90
            dcopy(NRHS, BX( I, 1 ), LDBX, B( PERM( I ), 1 ), LDB );
         } // 90

         // Step (4R): apply back the Givens rotations performed.

         DO 100 I = GIVPTR, 1, -1;
            drot(NRHS, B( GIVCOL( I, 2 ), 1 ), LDB, B( GIVCOL( I, 1 ), 1 ), LDB, GIVNUM( I, 2 ), -GIVNUM( I, 1 ) );
         } // 100
      }

      return;

      // End of DLALS0

      }
