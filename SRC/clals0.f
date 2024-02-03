      SUBROUTINE CLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B, LDB, BX, LDBX, PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, POLES, DIFL, DIFR, Z, K, C, S, RWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                GIVPTR, ICOMPQ, INFO, K, LDB, LDBX, LDGCOL, LDGNUM, NL, NR, NRHS, SQRE;
      REAL               C, S
      // ..
      // .. Array Arguments ..
      int                GIVCOL( LDGCOL, * ), PERM( * );
      REAL               DIFL( * ), DIFR( LDGNUM, * ), GIVNUM( LDGNUM, * ), POLES( LDGNUM, * ), RWORK( * ), Z( * )
      COMPLEX            B( LDB, * ), BX( LDBX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO, NEGONE
      const              ONE = 1.0E0, ZERO = 0.0E0, NEGONE = -1.0E0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, JCOL, JROW, M, N, NLP1;
      REAL               DIFLJ, DIFRJ, DJ, DSIGJ, DSIGJP, TEMP
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CLACPY, CLASCL, CSROT, CSSCAL, SGEMV, XERBLA
      // ..
      // .. External Functions ..
      REAL               SLAMC3, SNRM2
      // EXTERNAL SLAMC3, SNRM2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, CMPLX, MAX, REAL
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
         xerbla('CLALS0', -INFO );
         RETURN
      }

      M = N + SQRE
      NLP1 = NL + 1

      if ( ICOMPQ.EQ.0 ) {

         // Apply back orthogonal transformations from the left.

         // Step (1L): apply back the Givens rotations performed.

         for (I = 1; I <= GIVPTR; I++) { // 10
            csrot(NRHS, B( GIVCOL( I, 2 ), 1 ), LDB, B( GIVCOL( I, 1 ), 1 ), LDB, GIVNUM( I, 2 ), GIVNUM( I, 1 ) );
   10    CONTINUE

         // Step (2L): permute rows of B.

         ccopy(NRHS, B( NLP1, 1 ), LDB, BX( 1, 1 ), LDBX );
         for (I = 2; I <= N; I++) { // 20
            ccopy(NRHS, B( PERM( I ), 1 ), LDB, BX( I, 1 ), LDBX );
   20    CONTINUE

         // Step (3L): apply the inverse of the left singular vector
         // matrix to BX.

         if ( K.EQ.1 ) {
            ccopy(NRHS, BX, LDBX, B, LDB );
            if ( Z( 1 ).LT.ZERO ) {
               csscal(NRHS, NEGONE, B, LDB );
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
               if ( ( Z( J ).EQ.ZERO ) .OR. ( POLES( J, 2 ).EQ.ZERO ) ) {
                  RWORK( J ) = ZERO
               } else {
                  RWORK( J ) = -POLES( J, 2 )*Z( J ) / DIFLJ / ( POLES( J, 2 )+DJ )
               }
               DO 30 I = 1, J - 1
                  if ( ( Z( I ).EQ.ZERO ) .OR. ( POLES( I, 2 ).EQ.ZERO ) ) {
                     RWORK( I ) = ZERO
                  } else {

                     // Use calls to the subroutine SLAMC3 to enforce the
                     // parentheses (x+y)+z. The goal is to prevent
                     // optimizing compilers from doing x+(y+z).

                     RWORK( I ) = POLES( I, 2 )*Z( I ) / ( SLAMC3( POLES( I, 2 ), DSIGJ )- DIFLJ ) / ( POLES( I, 2 )+DJ )
                  }
   30          CONTINUE
               DO 40 I = J + 1, K
                  if ( ( Z( I ).EQ.ZERO ) .OR. ( POLES( I, 2 ).EQ.ZERO ) ) {
                     RWORK( I ) = ZERO
                  } else {
                     RWORK( I ) = POLES( I, 2 )*Z( I ) / ( SLAMC3( POLES( I, 2 ), DSIGJP )+ DIFRJ ) / ( POLES( I, 2 )+DJ )
                  }
   40          CONTINUE
               RWORK( 1 ) = NEGONE
               TEMP = SNRM2( K, RWORK, 1 )

               // Since B and BX are complex, the following call to SGEMV
               // is performed in two steps (real and imaginary parts).

               // CALL SGEMV( 'T', K, NRHS, ONE, BX, LDBX, WORK, 1, ZERO,
*    $                     B( J, 1 ), LDB )

               I = K + NRHS*2
               for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 60
                  for (JROW = 1; JROW <= K; JROW++) { // 50
                     I = I + 1
                     RWORK( I ) = REAL( BX( JROW, JCOL ) )
   50             CONTINUE
   60          CONTINUE
               sgemv('T', K, NRHS, ONE, RWORK( 1+K+NRHS*2 ), K, RWORK( 1 ), 1, ZERO, RWORK( 1+K ), 1 );
               I = K + NRHS*2
               for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 80
                  for (JROW = 1; JROW <= K; JROW++) { // 70
                     I = I + 1
                     RWORK( I ) = AIMAG( BX( JROW, JCOL ) )
   70             CONTINUE
   80          CONTINUE
               sgemv('T', K, NRHS, ONE, RWORK( 1+K+NRHS*2 ), K, RWORK( 1 ), 1, ZERO, RWORK( 1+K+NRHS ), 1 );
               for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 90
                  B( J, JCOL ) = CMPLX( RWORK( JCOL+K ), RWORK( JCOL+K+NRHS ) )
   90          CONTINUE
               clascl('G', 0, 0, TEMP, ONE, 1, NRHS, B( J, 1 ), LDB, INFO );
  100       CONTINUE
         }

         // Move the deflated rows of BX to B also.

         IF( K.LT.MAX( M, N ) ) CALL CLACPY( 'A', N-K, NRHS, BX( K+1, 1 ), LDBX, B( K+1, 1 ), LDB )
      } else {

         // Apply back the right orthogonal transformations.

         // Step (1R): apply back the new right singular vector matrix
         // to B.

         if ( K.EQ.1 ) {
            ccopy(NRHS, B, LDB, BX, LDBX );
         } else {
            for (J = 1; J <= K; J++) { // 180
               DSIGJ = POLES( J, 2 )
               if ( Z( J ).EQ.ZERO ) {
                  RWORK( J ) = ZERO
               } else {
                  RWORK( J ) = -Z( J ) / DIFL( J ) / ( DSIGJ+POLES( J, 1 ) ) / DIFR( J, 2 )
               }
               DO 110 I = 1, J - 1
                  if ( Z( J ).EQ.ZERO ) {
                     RWORK( I ) = ZERO
                  } else {

                     // Use calls to the subroutine SLAMC3 to enforce the
                     // parentheses (x+y)+z. The goal is to prevent optimizing
                     // compilers from doing x+(y+z).

                     RWORK( I ) = Z( J ) / ( SLAMC3( DSIGJ, -POLES( I+1, 2 ) )-DIFR( I, 1 ) ) / ( DSIGJ+POLES( I, 1 ) ) / DIFR( I, 2 )
                  }
  110          CONTINUE
               DO 120 I = J + 1, K
                  if ( Z( J ).EQ.ZERO ) {
                     RWORK( I ) = ZERO
                  } else {
                     RWORK( I ) = Z( J ) / ( SLAMC3( DSIGJ, -POLES( I, 2 ) )-DIFL( I ) ) / ( DSIGJ+POLES( I, 1 ) ) / DIFR( I, 2 )
                  }
  120          CONTINUE

               // Since B and BX are complex, the following call to SGEMV
               // is performed in two steps (real and imaginary parts).

               // CALL SGEMV( 'T', K, NRHS, ONE, B, LDB, WORK, 1, ZERO,
*    $                     BX( J, 1 ), LDBX )

               I = K + NRHS*2
               for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 140
                  for (JROW = 1; JROW <= K; JROW++) { // 130
                     I = I + 1
                     RWORK( I ) = REAL( B( JROW, JCOL ) )
  130             CONTINUE
  140          CONTINUE
               sgemv('T', K, NRHS, ONE, RWORK( 1+K+NRHS*2 ), K, RWORK( 1 ), 1, ZERO, RWORK( 1+K ), 1 );
               I = K + NRHS*2
               for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 160
                  for (JROW = 1; JROW <= K; JROW++) { // 150
                     I = I + 1
                     RWORK( I ) = AIMAG( B( JROW, JCOL ) )
  150             CONTINUE
  160          CONTINUE
               sgemv('T', K, NRHS, ONE, RWORK( 1+K+NRHS*2 ), K, RWORK( 1 ), 1, ZERO, RWORK( 1+K+NRHS ), 1 );
               for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 170
                  BX( J, JCOL ) = CMPLX( RWORK( JCOL+K ), RWORK( JCOL+K+NRHS ) )
  170          CONTINUE
  180       CONTINUE
         }

         // Step (2R): if SQRE = 1, apply back the rotation that is
         // related to the right null space of the subproblem.

         if ( SQRE.EQ.1 ) {
            ccopy(NRHS, B( M, 1 ), LDB, BX( M, 1 ), LDBX );
            csrot(NRHS, BX( 1, 1 ), LDBX, BX( M, 1 ), LDBX, C, S );
         }
         IF( K.LT.MAX( M, N ) ) CALL CLACPY( 'A', N-K, NRHS, B( K+1, 1 ), LDB, BX( K+1, 1 ), LDBX )

         // Step (3R): permute rows of B.

         ccopy(NRHS, BX( 1, 1 ), LDBX, B( NLP1, 1 ), LDB );
         if ( SQRE.EQ.1 ) {
            ccopy(NRHS, BX( M, 1 ), LDBX, B( M, 1 ), LDB );
         }
         for (I = 2; I <= N; I++) { // 190
            ccopy(NRHS, BX( I, 1 ), LDBX, B( PERM( I ), 1 ), LDB );
  190    CONTINUE

         // Step (4R): apply back the Givens rotations performed.

         DO 200 I = GIVPTR, 1, -1
            csrot(NRHS, B( GIVCOL( I, 2 ), 1 ), LDB, B( GIVCOL( I, 1 ), 1 ), LDB, GIVNUM( I, 2 ), -GIVNUM( I, 1 ) );
  200    CONTINUE
      }

      RETURN

      // End of CLALS0

      }
