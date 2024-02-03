      SUBROUTINE DLALSA( ICOMPQ, SMLSIZ, N, NRHS, B, LDB, BX, LDBX, U, LDU, VT, K, DIFL, DIFR, Z, POLES, GIVPTR, GIVCOL, LDGCOL, PERM, GIVNUM, C, S, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                ICOMPQ, INFO, LDB, LDBX, LDGCOL, LDU, N, NRHS, SMLSIZ;
      // ..
      // .. Array Arguments ..
      int                GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ), K( * ), PERM( LDGCOL, * );
      double             B( LDB, * ), BX( LDBX, * ), C( * ), DIFL( LDU, * ), DIFR( LDU, * ), GIVNUM( LDU, * ), POLES( LDU, * ), S( * ), U( LDU, * ), VT( LDU, * ), WORK( * ), Z( LDU, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      int                I, I1, IC, IM1, INODE, J, LF, LL, LVL, LVL2, ND, NDB1, NDIML, NDIMR, NL, NLF, NLP1, NLVL, NR, NRF, NRP1, SQRE;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGEMM, DLALS0, DLASDT, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0

      if ( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) {
         INFO = -1
      } else if ( SMLSIZ.LT.3 ) {
         INFO = -2
      } else if ( N.LT.SMLSIZ ) {
         INFO = -3
      } else if ( NRHS.LT.1 ) {
         INFO = -4
      } else if ( LDB.LT.N ) {
         INFO = -6
      } else if ( LDBX.LT.N ) {
         INFO = -8
      } else if ( LDU.LT.N ) {
         INFO = -10
      } else if ( LDGCOL.LT.N ) {
         INFO = -19
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'DLALSA', -INFO )
         RETURN
      }

      // Book-keeping and  setting up the computation tree.

      INODE = 1
      NDIML = INODE + N
      NDIMR = NDIML + N

      CALL DLASDT( N, NLVL, ND, IWORK( INODE ), IWORK( NDIML ), IWORK( NDIMR ), SMLSIZ )

      // The following code applies back the left singular vector factors.
      // For applying back the right singular vector factors, go to 50.

      if ( ICOMPQ.EQ.1 ) {
         GO TO 50
      }

      // The nodes on the bottom level of the tree were solved
      // by DLASDQ. The corresponding left and right singular vector
      // matrices are in explicit form. First apply back the left
      // singular vector matrices.

      NDB1 = ( ND+1 ) / 2
      DO 10 I = NDB1, ND

         // IC : center row of each node
         // NL : number of rows of left  subproblem
         // NR : number of rows of right subproblem
         // NLF: starting row of the left   subproblem
         // NRF: starting row of the right  subproblem

         I1 = I - 1
         IC = IWORK( INODE+I1 )
         NL = IWORK( NDIML+I1 )
         NR = IWORK( NDIMR+I1 )
         NLF = IC - NL
         NRF = IC + 1
         CALL DGEMM( 'T', 'N', NL, NRHS, NL, ONE, U( NLF, 1 ), LDU, B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX )          CALL DGEMM( 'T', 'N', NR, NRHS, NR, ONE, U( NRF, 1 ), LDU, B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX )
   10 CONTINUE

      // Next copy the rows of B that correspond to unchanged rows
      // in the bidiagonal matrix to BX.

      DO 20 I = 1, ND
         IC = IWORK( INODE+I-1 )
         CALL DCOPY( NRHS, B( IC, 1 ), LDB, BX( IC, 1 ), LDBX )
   20 CONTINUE

      // Finally go through the left singular vector matrices of all
     t // he other subproblems bottom-up on the tree.

      J = 2**NLVL
      SQRE = 0

      DO 40 LVL = NLVL, 1, -1
         LVL2 = 2*LVL - 1

         // find the first node LF and last node LL on
        t // he current level LVL

         if ( LVL.EQ.1 ) {
            LF = 1
            LL = 1
         } else {
            LF = 2**( LVL-1 )
            LL = 2*LF - 1
         }
         DO 30 I = LF, LL
            IM1 = I - 1
            IC = IWORK( INODE+IM1 )
            NL = IWORK( NDIML+IM1 )
            NR = IWORK( NDIMR+IM1 )
            NLF = IC - NL
            NRF = IC + 1
            J = J - 1
            CALL DLALS0( ICOMPQ, NL, NR, SQRE, NRHS, BX( NLF, 1 ), LDBX, B( NLF, 1 ), LDB, PERM( NLF, LVL ), GIVPTR( J ), GIVCOL( NLF, LVL2 ), LDGCOL, GIVNUM( NLF, LVL2 ), LDU, POLES( NLF, LVL2 ), DIFL( NLF, LVL ), DIFR( NLF, LVL2 ), Z( NLF, LVL ), K( J ), C( J ), S( J ), WORK, INFO )
   30    CONTINUE
   40 CONTINUE
      GO TO 90

      // ICOMPQ = 1: applying back the right singular vector factors.

   50 CONTINUE

      // First now go through the right singular vector matrices of all
     t // he tree nodes top-down.

      J = 0
      DO 70 LVL = 1, NLVL
         LVL2 = 2*LVL - 1

         // Find the first node LF and last node LL on
        t // he current level LVL.

         if ( LVL.EQ.1 ) {
            LF = 1
            LL = 1
         } else {
            LF = 2**( LVL-1 )
            LL = 2*LF - 1
         }
         DO 60 I = LL, LF, -1
            IM1 = I - 1
            IC = IWORK( INODE+IM1 )
            NL = IWORK( NDIML+IM1 )
            NR = IWORK( NDIMR+IM1 )
            NLF = IC - NL
            NRF = IC + 1
            if ( I.EQ.LL ) {
               SQRE = 0
            } else {
               SQRE = 1
            }
            J = J + 1
            CALL DLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B( NLF, 1 ), LDB, BX( NLF, 1 ), LDBX, PERM( NLF, LVL ), GIVPTR( J ), GIVCOL( NLF, LVL2 ), LDGCOL, GIVNUM( NLF, LVL2 ), LDU, POLES( NLF, LVL2 ), DIFL( NLF, LVL ), DIFR( NLF, LVL2 ), Z( NLF, LVL ), K( J ), C( J ), S( J ), WORK, INFO )
   60    CONTINUE
   70 CONTINUE

      // The nodes on the bottom level of the tree were solved
      // by DLASDQ. The corresponding right singular vector
      // matrices are in explicit form. Apply them back.

      NDB1 = ( ND+1 ) / 2
      DO 80 I = NDB1, ND
         I1 = I - 1
         IC = IWORK( INODE+I1 )
         NL = IWORK( NDIML+I1 )
         NR = IWORK( NDIMR+I1 )
         NLP1 = NL + 1
         if ( I.EQ.ND ) {
            NRP1 = NR
         } else {
            NRP1 = NR + 1
         }
         NLF = IC - NL
         NRF = IC + 1
         CALL DGEMM( 'T', 'N', NLP1, NRHS, NLP1, ONE, VT( NLF, 1 ), LDU, B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX )          CALL DGEMM( 'T', 'N', NRP1, NRHS, NRP1, ONE, VT( NRF, 1 ), LDU, B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX )
   80 CONTINUE

   90 CONTINUE

      RETURN

      // End of DLALSA

      }
