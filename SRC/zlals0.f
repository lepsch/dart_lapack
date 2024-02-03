      SUBROUTINE ZLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B, LDB, BX, LDBX, PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, POLES, DIFL, DIFR, Z, K, C, S, RWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                GIVPTR, ICOMPQ, INFO, K, LDB, LDBX, LDGCOL, LDGNUM, NL, NR, NRHS, SQRE;
      double             C, S;
*     ..
*     .. Array Arguments ..
      int                GIVCOL( LDGCOL, * ), PERM( * );
      double             DIFL( * ), DIFR( LDGNUM, * ), GIVNUM( LDGNUM, * ), POLES( LDGNUM, * ), RWORK( * ), Z( * );
      COMPLEX*16         B( LDB, * ), BX( LDBX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ONE, ZERO, NEGONE;
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0, NEGONE = -1.0D0 )
*     ..
*     .. Local Scalars ..
      int                I, J, JCOL, JROW, M, N, NLP1;
      double             DIFLJ, DIFRJ, DJ, DSIGJ, DSIGJP, TEMP;
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, XERBLA, ZCOPY, ZDROT, ZDSCAL, ZLACPY, ZLASCL
*     ..
*     .. External Functions ..
      double             DLAMC3, DNRM2;
      EXTERNAL           DLAMC3, DNRM2
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, DIMAG, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
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
         CALL XERBLA( 'ZLALS0', -INFO )
         RETURN
      END IF
*
      M = N + SQRE
      NLP1 = NL + 1
*
      IF( ICOMPQ.EQ.0 ) THEN
*
*        Apply back orthogonal transformations from the left.
*
*        Step (1L): apply back the Givens rotations performed.
*
         DO 10 I = 1, GIVPTR
            CALL ZDROT( NRHS, B( GIVCOL( I, 2 ), 1 ), LDB, B( GIVCOL( I, 1 ), 1 ), LDB, GIVNUM( I, 2 ), GIVNUM( I, 1 ) )
   10    CONTINUE
*
*        Step (2L): permute rows of B.
*
         CALL ZCOPY( NRHS, B( NLP1, 1 ), LDB, BX( 1, 1 ), LDBX )
         DO 20 I = 2, N
            CALL ZCOPY( NRHS, B( PERM( I ), 1 ), LDB, BX( I, 1 ), LDBX )
   20    CONTINUE
*
*        Step (3L): apply the inverse of the left singular vector
*        matrix to BX.
*
         IF( K.EQ.1 ) THEN
            CALL ZCOPY( NRHS, BX, LDBX, B, LDB )
            IF( Z( 1 ).LT.ZERO ) THEN
               CALL ZDSCAL( NRHS, NEGONE, B, LDB )
            END IF
         ELSE
            DO 100 J = 1, K
               DIFLJ = DIFL( J )
               DJ = POLES( J, 1 )
               DSIGJ = -POLES( J, 2 )
               IF( J.LT.K ) THEN
                  DIFRJ = -DIFR( J, 1 )
                  DSIGJP = -POLES( J+1, 2 )
               END IF
               IF( ( Z( J ).EQ.ZERO ) .OR. ( POLES( J, 2 ).EQ.ZERO ) ) THEN
                  RWORK( J ) = ZERO
               ELSE
                  RWORK( J ) = -POLES( J, 2 )*Z( J ) / DIFLJ / ( POLES( J, 2 )+DJ )
               END IF
               DO 30 I = 1, J - 1
                  IF( ( Z( I ).EQ.ZERO ) .OR. ( POLES( I, 2 ).EQ.ZERO ) ) THEN
                     RWORK( I ) = ZERO
                  ELSE
*
*                    Use calls to the subroutine DLAMC3 to enforce the
*                    parentheses (x+y)+z. The goal is to prevent
*                    optimizing compilers from doing x+(y+z).
*
                     RWORK( I ) = POLES( I, 2 )*Z( I ) / ( DLAMC3( POLES( I, 2 ), DSIGJ )- DIFLJ ) / ( POLES( I, 2 )+DJ )
                  END IF
   30          CONTINUE
               DO 40 I = J + 1, K
                  IF( ( Z( I ).EQ.ZERO ) .OR. ( POLES( I, 2 ).EQ.ZERO ) ) THEN
                     RWORK( I ) = ZERO
                  ELSE
                     RWORK( I ) = POLES( I, 2 )*Z( I ) / ( DLAMC3( POLES( I, 2 ), DSIGJP )+ DIFRJ ) / ( POLES( I, 2 )+DJ )
                  END IF
   40          CONTINUE
               RWORK( 1 ) = NEGONE
               TEMP = DNRM2( K, RWORK, 1 )
*
*              Since B and BX are complex, the following call to DGEMV
*              is performed in two steps (real and imaginary parts).
*
*              CALL DGEMV( 'T', K, NRHS, ONE, BX, LDBX, WORK, 1, ZERO,
*    $                     B( J, 1 ), LDB )
*
               I = K + NRHS*2
               DO 60 JCOL = 1, NRHS
                  DO 50 JROW = 1, K
                     I = I + 1
                     RWORK( I ) = DBLE( BX( JROW, JCOL ) )
   50             CONTINUE
   60          CONTINUE
               CALL DGEMV( 'T', K, NRHS, ONE, RWORK( 1+K+NRHS*2 ), K, RWORK( 1 ), 1, ZERO, RWORK( 1+K ), 1 )
               I = K + NRHS*2
               DO 80 JCOL = 1, NRHS
                  DO 70 JROW = 1, K
                     I = I + 1
                     RWORK( I ) = DIMAG( BX( JROW, JCOL ) )
   70             CONTINUE
   80          CONTINUE
               CALL DGEMV( 'T', K, NRHS, ONE, RWORK( 1+K+NRHS*2 ), K, RWORK( 1 ), 1, ZERO, RWORK( 1+K+NRHS ), 1 )
               DO 90 JCOL = 1, NRHS
                  B( J, JCOL ) = DCMPLX( RWORK( JCOL+K ), RWORK( JCOL+K+NRHS ) )
   90          CONTINUE
               CALL ZLASCL( 'G', 0, 0, TEMP, ONE, 1, NRHS, B( J, 1 ), LDB, INFO )
  100       CONTINUE
         END IF
*
*        Move the deflated rows of BX to B also.
*
         IF( K.LT.MAX( M, N ) ) CALL ZLACPY( 'A', N-K, NRHS, BX( K+1, 1 ), LDBX, B( K+1, 1 ), LDB )
      ELSE
*
*        Apply back the right orthogonal transformations.
*
*        Step (1R): apply back the new right singular vector matrix
*        to B.
*
         IF( K.EQ.1 ) THEN
            CALL ZCOPY( NRHS, B, LDB, BX, LDBX )
         ELSE
            DO 180 J = 1, K
               DSIGJ = POLES( J, 2 )
               IF( Z( J ).EQ.ZERO ) THEN
                  RWORK( J ) = ZERO
               ELSE
                  RWORK( J ) = -Z( J ) / DIFL( J ) / ( DSIGJ+POLES( J, 1 ) ) / DIFR( J, 2 )
               END IF
               DO 110 I = 1, J - 1
                  IF( Z( J ).EQ.ZERO ) THEN
                     RWORK( I ) = ZERO
                  ELSE
*
*                    Use calls to the subroutine DLAMC3 to enforce the
*                    parentheses (x+y)+z. The goal is to prevent
*                    optimizing compilers from doing x+(y+z).
*
                     RWORK( I ) = Z( J ) / ( DLAMC3( DSIGJ, -POLES( I+1, 2 ) )-DIFR( I, 1 ) ) / ( DSIGJ+POLES( I, 1 ) ) / DIFR( I, 2 )
                  END IF
  110          CONTINUE
               DO 120 I = J + 1, K
                  IF( Z( J ).EQ.ZERO ) THEN
                     RWORK( I ) = ZERO
                  ELSE
                     RWORK( I ) = Z( J ) / ( DLAMC3( DSIGJ, -POLES( I, 2 ) )-DIFL( I ) ) / ( DSIGJ+POLES( I, 1 ) ) / DIFR( I, 2 )
                  END IF
  120          CONTINUE
*
*              Since B and BX are complex, the following call to DGEMV
*              is performed in two steps (real and imaginary parts).
*
*              CALL DGEMV( 'T', K, NRHS, ONE, B, LDB, WORK, 1, ZERO,
*    $                     BX( J, 1 ), LDBX )
*
               I = K + NRHS*2
               DO 140 JCOL = 1, NRHS
                  DO 130 JROW = 1, K
                     I = I + 1
                     RWORK( I ) = DBLE( B( JROW, JCOL ) )
  130             CONTINUE
  140          CONTINUE
               CALL DGEMV( 'T', K, NRHS, ONE, RWORK( 1+K+NRHS*2 ), K, RWORK( 1 ), 1, ZERO, RWORK( 1+K ), 1 )
               I = K + NRHS*2
               DO 160 JCOL = 1, NRHS
                  DO 150 JROW = 1, K
                     I = I + 1
                     RWORK( I ) = DIMAG( B( JROW, JCOL ) )
  150             CONTINUE
  160          CONTINUE
               CALL DGEMV( 'T', K, NRHS, ONE, RWORK( 1+K+NRHS*2 ), K, RWORK( 1 ), 1, ZERO, RWORK( 1+K+NRHS ), 1 )
               DO 170 JCOL = 1, NRHS
                  BX( J, JCOL ) = DCMPLX( RWORK( JCOL+K ), RWORK( JCOL+K+NRHS ) )
  170          CONTINUE
  180       CONTINUE
         END IF
*
*        Step (2R): if SQRE = 1, apply back the rotation that is
*        related to the right null space of the subproblem.
*
         IF( SQRE.EQ.1 ) THEN
            CALL ZCOPY( NRHS, B( M, 1 ), LDB, BX( M, 1 ), LDBX )
            CALL ZDROT( NRHS, BX( 1, 1 ), LDBX, BX( M, 1 ), LDBX, C, S )
         END IF
         IF( K.LT.MAX( M, N ) ) CALL ZLACPY( 'A', N-K, NRHS, B( K+1, 1 ), LDB, BX( K+1, 1 ), LDBX )
*
*        Step (3R): permute rows of B.
*
         CALL ZCOPY( NRHS, BX( 1, 1 ), LDBX, B( NLP1, 1 ), LDB )
         IF( SQRE.EQ.1 ) THEN
            CALL ZCOPY( NRHS, BX( M, 1 ), LDBX, B( M, 1 ), LDB )
         END IF
         DO 190 I = 2, N
            CALL ZCOPY( NRHS, BX( I, 1 ), LDBX, B( PERM( I ), 1 ), LDB )
  190    CONTINUE
*
*        Step (4R): apply back the Givens rotations performed.
*
         DO 200 I = GIVPTR, 1, -1
            CALL ZDROT( NRHS, B( GIVCOL( I, 2 ), 1 ), LDB, B( GIVCOL( I, 1 ), 1 ), LDB, GIVNUM( I, 2 ), -GIVNUM( I, 1 ) )
  200    CONTINUE
      END IF
*
      RETURN
*
*     End of ZLALS0
*
      END
