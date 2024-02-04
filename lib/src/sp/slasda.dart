      void slasda(ICOMPQ, SMLSIZ, N, SQRE, D, E, U, LDU, VT, K, DIFL, DIFR, Z, POLES, GIVPTR, GIVCOL, LDGCOL, PERM, GIVNUM, C, S, WORK, IWORK, INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                ICOMPQ, INFO, LDGCOL, LDU, N, SMLSIZ, SQRE;
      // ..
      // .. Array Arguments ..
      int                GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ), K( * ), PERM( LDGCOL, * );
      REAL               C( * ), D( * ), DIFL( LDU, * ), DIFR( LDU, * ), E( * ), GIVNUM( LDU, * ), POLES( LDU, * ), S( * ), U( LDU, * ), VT( LDU, * ), WORK( * ), Z( LDU, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, I1, IC, IDXQ, IDXQI, IM1, INODE, ITEMP, IWK, J, LF, LL, LVL, LVL2, M, NCC, ND, NDB1, NDIML, NDIMR, NL, NLF, NLP1, NLVL, NR, NRF, NRP1, NRU, NWORK1, NWORK2, SMLSZP, SQREI, VF, VFI, VL, VLI;
      REAL               ALPHA, BETA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SLASD6, SLASDQ, SLASDT, SLASET, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;

      if ( ( ICOMPQ < 0 ) || ( ICOMPQ > 1 ) ) {
         INFO = -1;
      } else if ( SMLSIZ < 3 ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( ( SQRE < 0 ) || ( SQRE > 1 ) ) {
         INFO = -4;
      } else if ( LDU < ( N+SQRE ) ) {
         INFO = -8;
      } else if ( LDGCOL < N ) {
         INFO = -17;
      }
      if ( INFO != 0 ) {
         xerbla('SLASDA', -INFO );
         return;
      }

      M = N + SQRE;

      // If the input matrix is too small, call SLASDQ to find the SVD.

      if ( N <= SMLSIZ ) {
         if ( ICOMPQ == 0 ) {
            slasdq('U', SQRE, N, 0, 0, 0, D, E, VT, LDU, U, LDU, U, LDU, WORK, INFO );
         } else {
            slasdq('U', SQRE, N, M, N, 0, D, E, VT, LDU, U, LDU, U, LDU, WORK, INFO );
         }
         return;
      }

      // Book-keeping and  set up the computation tree.

      INODE = 1;
      NDIML = INODE + N;
      NDIMR = NDIML + N;
      IDXQ = NDIMR + N;
      IWK = IDXQ + N;

      NCC = 0;
      NRU = 0;

      SMLSZP = SMLSIZ + 1;
      VF = 1;
      VL = VF + M;
      NWORK1 = VL + M;
      NWORK2 = NWORK1 + SMLSZP*SMLSZP;

      slasdt(N, NLVL, ND, IWORK( INODE ), IWORK( NDIML ), IWORK( NDIMR ), SMLSIZ );

      // for the nodes on bottom level of the tree, solve
      // their subproblems by SLASDQ.

      NDB1 = ( ND+1 ) / 2;
      for (I = NDB1; I <= ND; I++) { // 30

         // IC : center row of each node
         // NL : number of rows of left  subproblem
         // NR : number of rows of right subproblem
         // NLF: starting row of the left   subproblem
         // NRF: starting row of the right  subproblem

         I1 = I - 1;
         IC = IWORK( INODE+I1 );
         NL = IWORK( NDIML+I1 );
         NLP1 = NL + 1;
         NR = IWORK( NDIMR+I1 );
         NLF = IC - NL;
         NRF = IC + 1;
         IDXQI = IDXQ + NLF - 2;
         VFI = VF + NLF - 1;
         VLI = VL + NLF - 1;
         SQREI = 1;
         if ( ICOMPQ == 0 ) {
            slaset('A', NLP1, NLP1, ZERO, ONE, WORK( NWORK1 ), SMLSZP );
            slasdq('U', SQREI, NL, NLP1, NRU, NCC, D( NLF ), E( NLF ), WORK( NWORK1 ), SMLSZP, WORK( NWORK2 ), NL, WORK( NWORK2 ), NL, WORK( NWORK2 ), INFO );
            ITEMP = NWORK1 + NL*SMLSZP;
            scopy(NLP1, WORK( NWORK1 ), 1, WORK( VFI ), 1 );
            scopy(NLP1, WORK( ITEMP ), 1, WORK( VLI ), 1 );
         } else {
            slaset('A', NL, NL, ZERO, ONE, U( NLF, 1 ), LDU );
            slaset('A', NLP1, NLP1, ZERO, ONE, VT( NLF, 1 ), LDU );
            slasdq('U', SQREI, NL, NLP1, NL, NCC, D( NLF ), E( NLF ), VT( NLF, 1 ), LDU, U( NLF, 1 ), LDU, U( NLF, 1 ), LDU, WORK( NWORK1 ), INFO );
            scopy(NLP1, VT( NLF, 1 ), 1, WORK( VFI ), 1 );
            scopy(NLP1, VT( NLF, NLP1 ), 1, WORK( VLI ), 1 );
         }
         if ( INFO != 0 ) {
            return;
         }
         for (J = 1; J <= NL; J++) { // 10
            IWORK[IDXQI+J] = J;
         } // 10
         if ( ( I == ND ) && ( SQRE == 0 ) ) {
            SQREI = 0;
         } else {
            SQREI = 1;
         }
         IDXQI = IDXQI + NLP1;
         VFI = VFI + NLP1;
         VLI = VLI + NLP1;
         NRP1 = NR + SQREI;
         if ( ICOMPQ == 0 ) {
            slaset('A', NRP1, NRP1, ZERO, ONE, WORK( NWORK1 ), SMLSZP );
            slasdq('U', SQREI, NR, NRP1, NRU, NCC, D( NRF ), E( NRF ), WORK( NWORK1 ), SMLSZP, WORK( NWORK2 ), NR, WORK( NWORK2 ), NR, WORK( NWORK2 ), INFO );
            ITEMP = NWORK1 + ( NRP1-1 )*SMLSZP;
            scopy(NRP1, WORK( NWORK1 ), 1, WORK( VFI ), 1 );
            scopy(NRP1, WORK( ITEMP ), 1, WORK( VLI ), 1 );
         } else {
            slaset('A', NR, NR, ZERO, ONE, U( NRF, 1 ), LDU );
            slaset('A', NRP1, NRP1, ZERO, ONE, VT( NRF, 1 ), LDU );
            slasdq('U', SQREI, NR, NRP1, NR, NCC, D( NRF ), E( NRF ), VT( NRF, 1 ), LDU, U( NRF, 1 ), LDU, U( NRF, 1 ), LDU, WORK( NWORK1 ), INFO );
            scopy(NRP1, VT( NRF, 1 ), 1, WORK( VFI ), 1 );
            scopy(NRP1, VT( NRF, NRP1 ), 1, WORK( VLI ), 1 );
         }
         if ( INFO != 0 ) {
            return;
         }
         for (J = 1; J <= NR; J++) { // 20
            IWORK[IDXQI+J] = J;
         } // 20
      } // 30

      // Now conquer each subproblem bottom-up.

      J = 2**NLVL;
      for (LVL = NLVL; LVL >= 1; LVL--) { // 50
         LVL2 = LVL*2 - 1;

         // Find the first node LF and last node LL on
         // the current level LVL.

         if ( LVL == 1 ) {
            LF = 1;
            LL = 1;
         } else {
            LF = 2**( LVL-1 );
            LL = 2*LF - 1;
         }
         for (I = LF; I <= LL; I++) { // 40
            IM1 = I - 1;
            IC = IWORK( INODE+IM1 );
            NL = IWORK( NDIML+IM1 );
            NR = IWORK( NDIMR+IM1 );
            NLF = IC - NL;
            NRF = IC + 1;
            if ( I == LL ) {
               SQREI = SQRE;
            } else {
               SQREI = 1;
            }
            VFI = VF + NLF - 1;
            VLI = VL + NLF - 1;
            IDXQI = IDXQ + NLF - 1;
            ALPHA = D( IC );
            BETA = E( IC );
            if ( ICOMPQ == 0 ) {
               slasd6(ICOMPQ, NL, NR, SQREI, D( NLF ), WORK( VFI ), WORK( VLI ), ALPHA, BETA, IWORK( IDXQI ), PERM, GIVPTR( 1 ), GIVCOL, LDGCOL, GIVNUM, LDU, POLES, DIFL, DIFR, Z, K( 1 ), C( 1 ), S( 1 ), WORK( NWORK1 ), IWORK( IWK ), INFO );
            } else {
               J = J - 1;
               slasd6(ICOMPQ, NL, NR, SQREI, D( NLF ), WORK( VFI ), WORK( VLI ), ALPHA, BETA, IWORK( IDXQI ), PERM( NLF, LVL ), GIVPTR( J ), GIVCOL( NLF, LVL2 ), LDGCOL, GIVNUM( NLF, LVL2 ), LDU, POLES( NLF, LVL2 ), DIFL( NLF, LVL ), DIFR( NLF, LVL2 ), Z( NLF, LVL ), K( J ), C( J ), S( J ), WORK( NWORK1 ), IWORK( IWK ), INFO );
            }
            if ( INFO != 0 ) {
               return;
            }
         } // 40
      } // 50

      return;
      }
