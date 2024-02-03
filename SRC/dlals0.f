      SUBROUTINE DLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B, LDB, BX, LDBX, PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, POLES, DIFL, DIFR, Z, K, C, S, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                GIVPTR, ICOMPQ, INFO, K, LDB, LDBX, LDGCOL, LDGNUM, NL, NR, NRHS, SQRE;
      double             C, S;
      // ..
      // .. Array Arguments ..
      int                GIVCOL( LDGCOL, * ), PERM( * );
      double             B( LDB, * ), BX( LDBX, * ), DIFL( * ), DIFR( LDGNUM, * ), GIVNUM( LDGNUM, * ), POLES( LDGNUM, * ), WORK( * ), Z( * );
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ONE, ZERO, NEGONE;
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0, NEGONE = -1.0D0 )
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
*
      // Test the input parameters.
*
      INFO = 0
      N = NL + NR + 1
*
      IF( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) THEN
         INFO = -1
      ELSE IF( NL.LT.1 ) THEN
         INFO = -2
      ELSE IF( NR.LT.1 ) THEN
         INFO = -3
      ELSE IF( ( SQRE.LT.0 ) .OR. ( SQRE.GT.1 ) ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.1 ) THEN
         INFO = -5
      ELSE IF( LDB.LT.N ) THEN
         INFO = -7
      ELSE IF( LDBX.LT.N ) THEN
         INFO = -9
      ELSE IF( GIVPTR.LT.0 ) THEN
         INFO = -11
      ELSE IF( LDGCOL.LT.N ) THEN
         INFO = -13
      ELSE IF( LDGNUM.LT.N ) THEN
         INFO = -15
      ELSE IF( K.LT.1 ) THEN
         INFO = -20
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLALS0', -INFO )
         RETURN
      END IF
*
      M = N + SQRE
      NLP1 = NL + 1
*
      IF( ICOMPQ.EQ.0 ) THEN
*
         // Apply back orthogonal transformations from the left.
*
         // Step (1L): apply back the Givens rotations performed.
*
         DO 10 I = 1, GIVPTR
            CALL DROT( NRHS, B( GIVCOL( I, 2 ), 1 ), LDB, B( GIVCOL( I, 1 ), 1 ), LDB, GIVNUM( I, 2 ), GIVNUM( I, 1 ) )
   10    CONTINUE
*
         // Step (2L): permute rows of B.
*
         CALL DCOPY( NRHS, B( NLP1, 1 ), LDB, BX( 1, 1 ), LDBX )
         DO 20 I = 2, N
            CALL DCOPY( NRHS, B( PERM( I ), 1 ), LDB, BX( I, 1 ), LDBX )
   20    CONTINUE
*
         // Step (3L): apply the inverse of the left singular vector
         // matrix to BX.
*
         IF( K.EQ.1 ) THEN
            CALL DCOPY( NRHS, BX, LDBX, B, LDB )
            IF( Z( 1 ).LT.ZERO ) THEN
               CALL DSCAL( NRHS, NEGONE, B, LDB )
            END IF
         ELSE
            DO 50 J = 1, K
               DIFLJ = DIFL( J )
               DJ = POLES( J, 1 )
               DSIGJ = -POLES( J, 2 )
               IF( J.LT.K ) THEN
                  DIFRJ = -DIFR( J, 1 )
                  DSIGJP = -POLES( J+1, 2 )
               END IF
               IF( ( Z( J ).EQ.ZERO ) .OR. ( POLES( J, 2 ).EQ.ZERO ) ) THEN
                  WORK( J ) = ZERO
               ELSE
                  WORK( J ) = -POLES( J, 2 )*Z( J ) / DIFLJ / ( POLES( J, 2 )+DJ )
               END IF
               DO 30 I = 1, J - 1
                  IF( ( Z( I ).EQ.ZERO ) .OR. ( POLES( I, 2 ).EQ.ZERO ) ) THEN
                     WORK( I ) = ZERO
                  ELSE
*
                     // Use calls to the subroutine DLAMC3 to enforce the
                     // parentheses (x+y)+z. The goal is to prevent
                     // optimizing compilers from doing x+(y+z).
*
                     WORK( I ) = POLES( I, 2 )*Z( I ) / ( DLAMC3( POLES( I, 2 ), DSIGJ )- DIFLJ ) / ( POLES( I, 2 )+DJ )
                  END IF
   30          CONTINUE
               DO 40 I = J + 1, K
                  IF( ( Z( I ).EQ.ZERO ) .OR. ( POLES( I, 2 ).EQ.ZERO ) ) THEN
                     WORK( I ) = ZERO
                  ELSE
                     WORK( I ) = POLES( I, 2 )*Z( I ) / ( DLAMC3( POLES( I, 2 ), DSIGJP )+ DIFRJ ) / ( POLES( I, 2 )+DJ )
                  END IF
   40          CONTINUE
               WORK( 1 ) = NEGONE
               TEMP = DNRM2( K, WORK, 1 )
               CALL DGEMV( 'T', K, NRHS, ONE, BX, LDBX, WORK, 1, ZERO, B( J, 1 ), LDB )                CALL DLASCL( 'G', 0, 0, TEMP, ONE, 1, NRHS, B( J, 1 ), LDB, INFO )
   50       CONTINUE
         END IF
*
         // Move the deflated rows of BX to B also.
*
         IF( K.LT.MAX( M, N ) ) CALL DLACPY( 'A', N-K, NRHS, BX( K+1, 1 ), LDBX, B( K+1, 1 ), LDB )
      ELSE
*
         // Apply back the right orthogonal transformations.
*
         // Step (1R): apply back the new right singular vector matrix
        t // o B.
*
         IF( K.EQ.1 ) THEN
            CALL DCOPY( NRHS, B, LDB, BX, LDBX )
         ELSE
            DO 80 J = 1, K
               DSIGJ = POLES( J, 2 )
               IF( Z( J ).EQ.ZERO ) THEN
                  WORK( J ) = ZERO
               ELSE
                  WORK( J ) = -Z( J ) / DIFL( J ) / ( DSIGJ+POLES( J, 1 ) ) / DIFR( J, 2 )
               END IF
               DO 60 I = 1, J - 1
                  IF( Z( J ).EQ.ZERO ) THEN
                     WORK( I ) = ZERO
                  ELSE
*
                     // Use calls to the subroutine DLAMC3 to enforce the
                     // parentheses (x+y)+z. The goal is to prevent
                     // optimizing compilers from doing x+(y+z).
*
                     WORK( I ) = Z( J ) / ( DLAMC3( DSIGJ, -POLES( I+1, 2 ) )-DIFR( I, 1 ) ) / ( DSIGJ+POLES( I, 1 ) ) / DIFR( I, 2 )
                  END IF
   60          CONTINUE
               DO 70 I = J + 1, K
                  IF( Z( J ).EQ.ZERO ) THEN
                     WORK( I ) = ZERO
                  ELSE
                     WORK( I ) = Z( J ) / ( DLAMC3( DSIGJ, -POLES( I, 2 ) )-DIFL( I ) ) / ( DSIGJ+POLES( I, 1 ) ) / DIFR( I, 2 )
                  END IF
   70          CONTINUE
               CALL DGEMV( 'T', K, NRHS, ONE, B, LDB, WORK, 1, ZERO, BX( J, 1 ), LDBX )
   80       CONTINUE
         END IF
*
         // Step (2R): if SQRE = 1, apply back the rotation that is
         // related to the right null space of the subproblem.
*
         IF( SQRE.EQ.1 ) THEN
            CALL DCOPY( NRHS, B( M, 1 ), LDB, BX( M, 1 ), LDBX )
            CALL DROT( NRHS, BX( 1, 1 ), LDBX, BX( M, 1 ), LDBX, C, S )
         END IF
         IF( K.LT.MAX( M, N ) ) CALL DLACPY( 'A', N-K, NRHS, B( K+1, 1 ), LDB, BX( K+1, 1 ), LDBX )
*
         // Step (3R): permute rows of B.
*
         CALL DCOPY( NRHS, BX( 1, 1 ), LDBX, B( NLP1, 1 ), LDB )
         IF( SQRE.EQ.1 ) THEN
            CALL DCOPY( NRHS, BX( M, 1 ), LDBX, B( M, 1 ), LDB )
         END IF
         DO 90 I = 2, N
            CALL DCOPY( NRHS, BX( I, 1 ), LDBX, B( PERM( I ), 1 ), LDB )
   90    CONTINUE
*
         // Step (4R): apply back the Givens rotations performed.
*
         DO 100 I = GIVPTR, 1, -1
            CALL DROT( NRHS, B( GIVCOL( I, 2 ), 1 ), LDB, B( GIVCOL( I, 1 ), 1 ), LDB, GIVNUM( I, 2 ), -GIVNUM( I, 1 ) )
  100    CONTINUE
      END IF
*
      RETURN
*
      // End of DLALS0
*
      END
