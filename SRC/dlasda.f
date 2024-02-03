      SUBROUTINE DLASDA( ICOMPQ, SMLSIZ, N, SQRE, D, E, U, LDU, VT, K, DIFL, DIFR, Z, POLES, GIVPTR, GIVCOL, LDGCOL, PERM, GIVNUM, C, S, WORK, IWORK, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                ICOMPQ, INFO, LDGCOL, LDU, N, SMLSIZ, SQRE;
      // ..
      // .. Array Arguments ..
      int                GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ), K( * ), PERM( LDGCOL, * )       double             C( * ), D( * ), DIFL( LDU, * ), DIFR( LDU, * ), E( * ), GIVNUM( LDU, * ), POLES( LDU, * ), S( * ), U( LDU, * ), VT( LDU, * ), WORK( * ), Z( LDU, * );;
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      // ..
      // .. Local Scalars ..
      int                I, I1, IC, IDXQ, IDXQI, IM1, INODE, ITEMP, IWK, J, LF, LL, LVL, LVL2, M, NCC, ND, NDB1, NDIML, NDIMR, NL, NLF, NLP1, NLVL, NR, NRF, NRP1, NRU, NWORK1, NWORK2, SMLSZP, SQREI, VF, VFI, VL, VLI;
      double             ALPHA, BETA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DLASD6, DLASDQ, DLASDT, DLASET, XERBLA
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters.
*
      INFO = 0
*
      IF( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) THEN
         INFO = -1
      ELSE IF( SMLSIZ.LT.3 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ( SQRE.LT.0 ) .OR. ( SQRE.GT.1 ) ) THEN
         INFO = -4
      ELSE IF( LDU.LT.( N+SQRE ) ) THEN
         INFO = -8
      ELSE IF( LDGCOL.LT.N ) THEN
         INFO = -17
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASDA', -INFO )
         RETURN
      END IF
*
      M = N + SQRE
*
      // If the input matrix is too small, call DLASDQ to find the SVD.
*
      IF( N.LE.SMLSIZ ) THEN
         IF( ICOMPQ.EQ.0 ) THEN
            CALL DLASDQ( 'U', SQRE, N, 0, 0, 0, D, E, VT, LDU, U, LDU, U, LDU, WORK, INFO )
         ELSE
            CALL DLASDQ( 'U', SQRE, N, M, N, 0, D, E, VT, LDU, U, LDU, U, LDU, WORK, INFO )
         END IF
         RETURN
      END IF
*
      // Book-keeping and  set up the computation tree.
*
      INODE = 1
      NDIML = INODE + N
      NDIMR = NDIML + N
      IDXQ = NDIMR + N
      IWK = IDXQ + N
*
      NCC = 0
      NRU = 0
*
      SMLSZP = SMLSIZ + 1
      VF = 1
      VL = VF + M
      NWORK1 = VL + M
      NWORK2 = NWORK1 + SMLSZP*SMLSZP
*
      CALL DLASDT( N, NLVL, ND, IWORK( INODE ), IWORK( NDIML ), IWORK( NDIMR ), SMLSIZ )
*
      // for the nodes on bottom level of the tree, solve
     t // heir subproblems by DLASDQ.
*
      NDB1 = ( ND+1 ) / 2
      DO 30 I = NDB1, ND
*
         // IC : center row of each node
         // NL : number of rows of left  subproblem
         // NR : number of rows of right subproblem
         // NLF: starting row of the left   subproblem
         // NRF: starting row of the right  subproblem
*
         I1 = I - 1
         IC = IWORK( INODE+I1 )
         NL = IWORK( NDIML+I1 )
         NLP1 = NL + 1
         NR = IWORK( NDIMR+I1 )
         NLF = IC - NL
         NRF = IC + 1
         IDXQI = IDXQ + NLF - 2
         VFI = VF + NLF - 1
         VLI = VL + NLF - 1
         SQREI = 1
         IF( ICOMPQ.EQ.0 ) THEN
            CALL DLASET( 'A', NLP1, NLP1, ZERO, ONE, WORK( NWORK1 ), SMLSZP )             CALL DLASDQ( 'U', SQREI, NL, NLP1, NRU, NCC, D( NLF ), E( NLF ), WORK( NWORK1 ), SMLSZP, WORK( NWORK2 ), NL, WORK( NWORK2 ), NL, WORK( NWORK2 ), INFO )
            ITEMP = NWORK1 + NL*SMLSZP
            CALL DCOPY( NLP1, WORK( NWORK1 ), 1, WORK( VFI ), 1 )
            CALL DCOPY( NLP1, WORK( ITEMP ), 1, WORK( VLI ), 1 )
         ELSE
            CALL DLASET( 'A', NL, NL, ZERO, ONE, U( NLF, 1 ), LDU )
            CALL DLASET( 'A', NLP1, NLP1, ZERO, ONE, VT( NLF, 1 ), LDU )
            CALL DLASDQ( 'U', SQREI, NL, NLP1, NL, NCC, D( NLF ), E( NLF ), VT( NLF, 1 ), LDU, U( NLF, 1 ), LDU, U( NLF, 1 ), LDU, WORK( NWORK1 ), INFO )
            CALL DCOPY( NLP1, VT( NLF, 1 ), 1, WORK( VFI ), 1 )
            CALL DCOPY( NLP1, VT( NLF, NLP1 ), 1, WORK( VLI ), 1 )
         END IF
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
         DO 10 J = 1, NL
            IWORK( IDXQI+J ) = J
   10    CONTINUE
         IF( ( I.EQ.ND ) .AND. ( SQRE.EQ.0 ) ) THEN
            SQREI = 0
         ELSE
            SQREI = 1
         END IF
         IDXQI = IDXQI + NLP1
         VFI = VFI + NLP1
         VLI = VLI + NLP1
         NRP1 = NR + SQREI
         IF( ICOMPQ.EQ.0 ) THEN
            CALL DLASET( 'A', NRP1, NRP1, ZERO, ONE, WORK( NWORK1 ), SMLSZP )             CALL DLASDQ( 'U', SQREI, NR, NRP1, NRU, NCC, D( NRF ), E( NRF ), WORK( NWORK1 ), SMLSZP, WORK( NWORK2 ), NR, WORK( NWORK2 ), NR, WORK( NWORK2 ), INFO )
            ITEMP = NWORK1 + ( NRP1-1 )*SMLSZP
            CALL DCOPY( NRP1, WORK( NWORK1 ), 1, WORK( VFI ), 1 )
            CALL DCOPY( NRP1, WORK( ITEMP ), 1, WORK( VLI ), 1 )
         ELSE
            CALL DLASET( 'A', NR, NR, ZERO, ONE, U( NRF, 1 ), LDU )
            CALL DLASET( 'A', NRP1, NRP1, ZERO, ONE, VT( NRF, 1 ), LDU )
            CALL DLASDQ( 'U', SQREI, NR, NRP1, NR, NCC, D( NRF ), E( NRF ), VT( NRF, 1 ), LDU, U( NRF, 1 ), LDU, U( NRF, 1 ), LDU, WORK( NWORK1 ), INFO )
            CALL DCOPY( NRP1, VT( NRF, 1 ), 1, WORK( VFI ), 1 )
            CALL DCOPY( NRP1, VT( NRF, NRP1 ), 1, WORK( VLI ), 1 )
         END IF
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
         DO 20 J = 1, NR
            IWORK( IDXQI+J ) = J
   20    CONTINUE
   30 CONTINUE
*
      // Now conquer each subproblem bottom-up.
*
      J = 2**NLVL
      DO 50 LVL = NLVL, 1, -1
         LVL2 = LVL*2 - 1
*
         // Find the first node LF and last node LL on
        t // he current level LVL.
*
         IF( LVL.EQ.1 ) THEN
            LF = 1
            LL = 1
         ELSE
            LF = 2**( LVL-1 )
            LL = 2*LF - 1
         END IF
         DO 40 I = LF, LL
            IM1 = I - 1
            IC = IWORK( INODE+IM1 )
            NL = IWORK( NDIML+IM1 )
            NR = IWORK( NDIMR+IM1 )
            NLF = IC - NL
            NRF = IC + 1
            IF( I.EQ.LL ) THEN
               SQREI = SQRE
            ELSE
               SQREI = 1
            END IF
            VFI = VF + NLF - 1
            VLI = VL + NLF - 1
            IDXQI = IDXQ + NLF - 1
            ALPHA = D( IC )
            BETA = E( IC )
            IF( ICOMPQ.EQ.0 ) THEN
               CALL DLASD6( ICOMPQ, NL, NR, SQREI, D( NLF ), WORK( VFI ), WORK( VLI ), ALPHA, BETA, IWORK( IDXQI ), PERM, GIVPTR( 1 ), GIVCOL, LDGCOL, GIVNUM, LDU, POLES, DIFL, DIFR, Z, K( 1 ), C( 1 ), S( 1 ), WORK( NWORK1 ), IWORK( IWK ), INFO )
            ELSE
               J = J - 1
               CALL DLASD6( ICOMPQ, NL, NR, SQREI, D( NLF ), WORK( VFI ), WORK( VLI ), ALPHA, BETA, IWORK( IDXQI ), PERM( NLF, LVL ), GIVPTR( J ), GIVCOL( NLF, LVL2 ), LDGCOL, GIVNUM( NLF, LVL2 ), LDU, POLES( NLF, LVL2 ), DIFL( NLF, LVL ), DIFR( NLF, LVL2 ), Z( NLF, LVL ), K( J ), C( J ), S( J ), WORK( NWORK1 ), IWORK( IWK ), INFO )
            END IF
            IF( INFO.NE.0 ) THEN
               RETURN
            END IF
   40    CONTINUE
   50 CONTINUE
*
      RETURN
*
      // End of DLASDA
*
      END
