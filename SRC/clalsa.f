      SUBROUTINE CLALSA( ICOMPQ, SMLSIZ, N, NRHS, B, LDB, BX, LDBX, U, LDU, VT, K, DIFL, DIFR, Z, POLES, GIVPTR, GIVCOL, LDGCOL, PERM, GIVNUM, C, S, RWORK, IWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                ICOMPQ, INFO, LDB, LDBX, LDGCOL, LDU, N, NRHS, SMLSIZ;
*     ..
*     .. Array Arguments ..
      int                GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ), K( * ), PERM( LDGCOL, * )       REAL               C( * ), DIFL( LDU, * ), DIFR( LDU, * ), GIVNUM( LDU, * ), POLES( LDU, * ), RWORK( * ), S( * ), U( LDU, * ), VT( LDU, * ), Z( LDU, * );
      COMPLEX            B( LDB, * ), BX( LDBX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      int                I, I1, IC, IM1, INODE, J, JCOL, JIMAG, JREAL, JROW, LF, LL, LVL, LVL2, ND, NDB1, NDIML, NDIMR, NL, NLF, NLP1, NLVL, NR, NRF, NRP1, SQRE;
*     ..
*     .. External Subroutines ..
      EXTERNAL           CCOPY, CLALS0, SGEMM, SLASDT, XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC AIMAG, CMPLX, REAL
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) THEN
         INFO = -1
      ELSE IF( SMLSIZ.LT.3 ) THEN
         INFO = -2
      ELSE IF( N.LT.SMLSIZ ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.1 ) THEN
         INFO = -4
      ELSE IF( LDB.LT.N ) THEN
         INFO = -6
      ELSE IF( LDBX.LT.N ) THEN
         INFO = -8
      ELSE IF( LDU.LT.N ) THEN
         INFO = -10
      ELSE IF( LDGCOL.LT.N ) THEN
         INFO = -19
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CLALSA', -INFO )
         RETURN
      END IF
*
*     Book-keeping and  setting up the computation tree.
*
      INODE = 1
      NDIML = INODE + N
      NDIMR = NDIML + N
*
      CALL SLASDT( N, NLVL, ND, IWORK( INODE ), IWORK( NDIML ), IWORK( NDIMR ), SMLSIZ )
*
*     The following code applies back the left singular vector factors.
*     For applying back the right singular vector factors, go to 170.
*
      IF( ICOMPQ.EQ.1 ) THEN
         GO TO 170
      END IF
*
*     The nodes on the bottom level of the tree were solved
*     by SLASDQ. The corresponding left and right singular vector
*     matrices are in explicit form. First apply back the left
*     singular vector matrices.
*
      NDB1 = ( ND+1 ) / 2
      DO 130 I = NDB1, ND
*
*        IC : center row of each node
*        NL : number of rows of left  subproblem
*        NR : number of rows of right subproblem
*        NLF: starting row of the left   subproblem
*        NRF: starting row of the right  subproblem
*
         I1 = I - 1
         IC = IWORK( INODE+I1 )
         NL = IWORK( NDIML+I1 )
         NR = IWORK( NDIMR+I1 )
         NLF = IC - NL
         NRF = IC + 1
*
*        Since B and BX are complex, the following call to SGEMM
*        is performed in two steps (real and imaginary parts).
*
*        CALL SGEMM( 'T', 'N', NL, NRHS, NL, ONE, U( NLF, 1 ), LDU,
*     $               B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX )
*
         J = NL*NRHS*2
         DO 20 JCOL = 1, NRHS
            DO 10 JROW = NLF, NLF + NL - 1
               J = J + 1
               RWORK( J ) = REAL( B( JROW, JCOL ) )
   10       CONTINUE
   20    CONTINUE
         CALL SGEMM( 'T', 'N', NL, NRHS, NL, ONE, U( NLF, 1 ), LDU, RWORK( 1+NL*NRHS*2 ), NL, ZERO, RWORK( 1 ), NL )
         J = NL*NRHS*2
         DO 40 JCOL = 1, NRHS
            DO 30 JROW = NLF, NLF + NL - 1
               J = J + 1
               RWORK( J ) = AIMAG( B( JROW, JCOL ) )
   30       CONTINUE
   40    CONTINUE
         CALL SGEMM( 'T', 'N', NL, NRHS, NL, ONE, U( NLF, 1 ), LDU, RWORK( 1+NL*NRHS*2 ), NL, ZERO, RWORK( 1+NL*NRHS ), NL )
         JREAL = 0
         JIMAG = NL*NRHS
         DO 60 JCOL = 1, NRHS
            DO 50 JROW = NLF, NLF + NL - 1
               JREAL = JREAL + 1
               JIMAG = JIMAG + 1
               BX( JROW, JCOL ) = CMPLX( RWORK( JREAL ), RWORK( JIMAG ) )
   50       CONTINUE
   60    CONTINUE
*
*        Since B and BX are complex, the following call to SGEMM
*        is performed in two steps (real and imaginary parts).
*
*        CALL SGEMM( 'T', 'N', NR, NRHS, NR, ONE, U( NRF, 1 ), LDU,
*    $               B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX )
*
         J = NR*NRHS*2
         DO 80 JCOL = 1, NRHS
            DO 70 JROW = NRF, NRF + NR - 1
               J = J + 1
               RWORK( J ) = REAL( B( JROW, JCOL ) )
   70       CONTINUE
   80    CONTINUE
         CALL SGEMM( 'T', 'N', NR, NRHS, NR, ONE, U( NRF, 1 ), LDU, RWORK( 1+NR*NRHS*2 ), NR, ZERO, RWORK( 1 ), NR )
         J = NR*NRHS*2
         DO 100 JCOL = 1, NRHS
            DO 90 JROW = NRF, NRF + NR - 1
               J = J + 1
               RWORK( J ) = AIMAG( B( JROW, JCOL ) )
   90       CONTINUE
  100    CONTINUE
         CALL SGEMM( 'T', 'N', NR, NRHS, NR, ONE, U( NRF, 1 ), LDU, RWORK( 1+NR*NRHS*2 ), NR, ZERO, RWORK( 1+NR*NRHS ), NR )
         JREAL = 0
         JIMAG = NR*NRHS
         DO 120 JCOL = 1, NRHS
            DO 110 JROW = NRF, NRF + NR - 1
               JREAL = JREAL + 1
               JIMAG = JIMAG + 1
               BX( JROW, JCOL ) = CMPLX( RWORK( JREAL ), RWORK( JIMAG ) )
  110       CONTINUE
  120    CONTINUE
*
  130 CONTINUE
*
*     Next copy the rows of B that correspond to unchanged rows
*     in the bidiagonal matrix to BX.
*
      DO 140 I = 1, ND
         IC = IWORK( INODE+I-1 )
         CALL CCOPY( NRHS, B( IC, 1 ), LDB, BX( IC, 1 ), LDBX )
  140 CONTINUE
*
*     Finally go through the left singular vector matrices of all
*     the other subproblems bottom-up on the tree.
*
      J = 2**NLVL
      SQRE = 0
*
      DO 160 LVL = NLVL, 1, -1
         LVL2 = 2*LVL - 1
*
*        find the first node LF and last node LL on
*        the current level LVL
*
         IF( LVL.EQ.1 ) THEN
            LF = 1
            LL = 1
         ELSE
            LF = 2**( LVL-1 )
            LL = 2*LF - 1
         END IF
         DO 150 I = LF, LL
            IM1 = I - 1
            IC = IWORK( INODE+IM1 )
            NL = IWORK( NDIML+IM1 )
            NR = IWORK( NDIMR+IM1 )
            NLF = IC - NL
            NRF = IC + 1
            J = J - 1
            CALL CLALS0( ICOMPQ, NL, NR, SQRE, NRHS, BX( NLF, 1 ), LDBX, B( NLF, 1 ), LDB, PERM( NLF, LVL ), GIVPTR( J ), GIVCOL( NLF, LVL2 ), LDGCOL, GIVNUM( NLF, LVL2 ), LDU, POLES( NLF, LVL2 ), DIFL( NLF, LVL ), DIFR( NLF, LVL2 ), Z( NLF, LVL ), K( J ), C( J ), S( J ), RWORK, INFO )
  150    CONTINUE
  160 CONTINUE
      GO TO 330
*
*     ICOMPQ = 1: applying back the right singular vector factors.
*
  170 CONTINUE
*
*     First now go through the right singular vector matrices of all
*     the tree nodes top-down.
*
      J = 0
      DO 190 LVL = 1, NLVL
         LVL2 = 2*LVL - 1
*
*        Find the first node LF and last node LL on
*        the current level LVL.
*
         IF( LVL.EQ.1 ) THEN
            LF = 1
            LL = 1
         ELSE
            LF = 2**( LVL-1 )
            LL = 2*LF - 1
         END IF
         DO 180 I = LL, LF, -1
            IM1 = I - 1
            IC = IWORK( INODE+IM1 )
            NL = IWORK( NDIML+IM1 )
            NR = IWORK( NDIMR+IM1 )
            NLF = IC - NL
            NRF = IC + 1
            IF( I.EQ.LL ) THEN
               SQRE = 0
            ELSE
               SQRE = 1
            END IF
            J = J + 1
            CALL CLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B( NLF, 1 ), LDB, BX( NLF, 1 ), LDBX, PERM( NLF, LVL ), GIVPTR( J ), GIVCOL( NLF, LVL2 ), LDGCOL, GIVNUM( NLF, LVL2 ), LDU, POLES( NLF, LVL2 ), DIFL( NLF, LVL ), DIFR( NLF, LVL2 ), Z( NLF, LVL ), K( J ), C( J ), S( J ), RWORK, INFO )
  180    CONTINUE
  190 CONTINUE
*
*     The nodes on the bottom level of the tree were solved
*     by SLASDQ. The corresponding right singular vector
*     matrices are in explicit form. Apply them back.
*
      NDB1 = ( ND+1 ) / 2
      DO 320 I = NDB1, ND
         I1 = I - 1
         IC = IWORK( INODE+I1 )
         NL = IWORK( NDIML+I1 )
         NR = IWORK( NDIMR+I1 )
         NLP1 = NL + 1
         IF( I.EQ.ND ) THEN
            NRP1 = NR
         ELSE
            NRP1 = NR + 1
         END IF
         NLF = IC - NL
         NRF = IC + 1
*
*        Since B and BX are complex, the following call to SGEMM is
*        performed in two steps (real and imaginary parts).
*
*        CALL SGEMM( 'T', 'N', NLP1, NRHS, NLP1, ONE, VT( NLF, 1 ), LDU,
*    $               B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX )
*
         J = NLP1*NRHS*2
         DO 210 JCOL = 1, NRHS
            DO 200 JROW = NLF, NLF + NLP1 - 1
               J = J + 1
               RWORK( J ) = REAL( B( JROW, JCOL ) )
  200       CONTINUE
  210    CONTINUE
         CALL SGEMM( 'T', 'N', NLP1, NRHS, NLP1, ONE, VT( NLF, 1 ), LDU, RWORK( 1+NLP1*NRHS*2 ), NLP1, ZERO, RWORK( 1 ), NLP1 )
         J = NLP1*NRHS*2
         DO 230 JCOL = 1, NRHS
            DO 220 JROW = NLF, NLF + NLP1 - 1
               J = J + 1
               RWORK( J ) = AIMAG( B( JROW, JCOL ) )
  220       CONTINUE
  230    CONTINUE
         CALL SGEMM( 'T', 'N', NLP1, NRHS, NLP1, ONE, VT( NLF, 1 ), LDU, RWORK( 1+NLP1*NRHS*2 ), NLP1, ZERO, RWORK( 1+NLP1*NRHS ), NLP1 )
         JREAL = 0
         JIMAG = NLP1*NRHS
         DO 250 JCOL = 1, NRHS
            DO 240 JROW = NLF, NLF + NLP1 - 1
               JREAL = JREAL + 1
               JIMAG = JIMAG + 1
               BX( JROW, JCOL ) = CMPLX( RWORK( JREAL ), RWORK( JIMAG ) )
  240       CONTINUE
  250    CONTINUE
*
*        Since B and BX are complex, the following call to SGEMM is
*        performed in two steps (real and imaginary parts).
*
*        CALL SGEMM( 'T', 'N', NRP1, NRHS, NRP1, ONE, VT( NRF, 1 ), LDU,
*    $               B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX )
*
         J = NRP1*NRHS*2
         DO 270 JCOL = 1, NRHS
            DO 260 JROW = NRF, NRF + NRP1 - 1
               J = J + 1
               RWORK( J ) = REAL( B( JROW, JCOL ) )
  260       CONTINUE
  270    CONTINUE
         CALL SGEMM( 'T', 'N', NRP1, NRHS, NRP1, ONE, VT( NRF, 1 ), LDU, RWORK( 1+NRP1*NRHS*2 ), NRP1, ZERO, RWORK( 1 ), NRP1 )
         J = NRP1*NRHS*2
         DO 290 JCOL = 1, NRHS
            DO 280 JROW = NRF, NRF + NRP1 - 1
               J = J + 1
               RWORK( J ) = AIMAG( B( JROW, JCOL ) )
  280       CONTINUE
  290    CONTINUE
         CALL SGEMM( 'T', 'N', NRP1, NRHS, NRP1, ONE, VT( NRF, 1 ), LDU, RWORK( 1+NRP1*NRHS*2 ), NRP1, ZERO, RWORK( 1+NRP1*NRHS ), NRP1 )
         JREAL = 0
         JIMAG = NRP1*NRHS
         DO 310 JCOL = 1, NRHS
            DO 300 JROW = NRF, NRF + NRP1 - 1
               JREAL = JREAL + 1
               JIMAG = JIMAG + 1
               BX( JROW, JCOL ) = CMPLX( RWORK( JREAL ), RWORK( JIMAG ) )
  300       CONTINUE
  310    CONTINUE
*
  320 CONTINUE
*
  330 CONTINUE
*
      RETURN
*
*     End of CLALSA
*
      END
