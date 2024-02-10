      void slasd0(N, SQRE, D, E, final Matrix<double> U, final int LDU, final Matrix<double> VT, final int LDVT, SMLSIZ, IWORK, WORK, Box<int> INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDU, LDVT, N, SMLSIZ, SQRE;
      int                IWORK( * );
      double               D( * ), E( * ), U( LDU, * ), VT( LDVT, * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, I1, IC, IDXQ, IDXQC, IM1, INODE, ITEMP, IWK, J, LF, LL, LVL, M, NCC, ND, NDB1, NDIML, NDIMR, NL, NLF, NLP1, NLVL, NR, NRF, NRP1, SQREI;
      double               ALPHA, BETA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASD1, SLASDQ, SLASDT, XERBLA

      // Test the input parameters.

      INFO = 0;

      if ( N < 0 ) {
         INFO = -1;
      } else if ( ( SQRE < 0 ) || ( SQRE > 1 ) ) {
         INFO = -2;
      }

      M = N + SQRE;

      if ( LDU < N ) {
         INFO = -6;
      } else if ( LDVT < M ) {
         INFO = -8;
      } else if ( SMLSIZ < 3 ) {
         INFO = -9;
      }
      if ( INFO != 0 ) {
         xerbla('SLASD0', -INFO );
         return;
      }

      // If the input matrix is too small, call SLASDQ to find the SVD.

      if ( N <= SMLSIZ ) {
         slasdq('U', SQRE, N, M, N, 0, D, E, VT, LDVT, U, LDU, U, LDU, WORK, INFO );
         return;
      }

      // Set up the computation tree.

      INODE = 1;
      NDIML = INODE + N;
      NDIMR = NDIML + N;
      IDXQ = NDIMR + N;
      IWK = IDXQ + N;
      slasdt(N, NLVL, ND, IWORK( INODE ), IWORK( NDIML ), IWORK( NDIMR ), SMLSIZ );

      // For the nodes on bottom level of the tree, solve
      // their subproblems by SLASDQ.

      NDB1 = ( ND+1 ) / 2;
      NCC = 0;
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
         NRP1 = NR + 1;
         NLF = IC - NL;
         NRF = IC + 1;
         SQREI = 1;
         slasdq('U', SQREI, NL, NLP1, NL, NCC, D( NLF ), E( NLF ), VT( NLF, NLF ), LDVT, U( NLF, NLF ), LDU, U( NLF, NLF ), LDU, WORK, INFO );
         if ( INFO != 0 ) {
            return;
         }
         ITEMP = IDXQ + NLF - 2;
         for (J = 1; J <= NL; J++) { // 10
            IWORK[ITEMP+J] = J;
         } // 10
         if ( I == ND ) {
            SQREI = SQRE;
         } else {
            SQREI = 1;
         }
         NRP1 = NR + SQREI;
         slasdq('U', SQREI, NR, NRP1, NR, NCC, D( NRF ), E( NRF ), VT( NRF, NRF ), LDVT, U( NRF, NRF ), LDU, U( NRF, NRF ), LDU, WORK, INFO );
         if ( INFO != 0 ) {
            return;
         }
         ITEMP = IDXQ + IC;
         for (J = 1; J <= NR; J++) { // 20
            IWORK[ITEMP+J-1] = J;
         } // 20
      } // 30

      // Now conquer each subproblem bottom-up.

      for (LVL = NLVL; LVL >= 1; LVL--) { // 50

         // Find the first node LF and last node LL on the
         // current level LVL.

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
            if ( ( SQRE == 0 ) && ( I == LL ) ) {
               SQREI = SQRE;
            } else {
               SQREI = 1;
            }
            IDXQC = IDXQ + NLF - 1;
            ALPHA = D( IC );
            BETA = E( IC );
            slasd1(NL, NR, SQREI, D( NLF ), ALPHA, BETA, U( NLF, NLF ), LDU, VT( NLF, NLF ), LDVT, IWORK( IDXQC ), IWORK( IWK ), WORK, INFO );

      // Report the possible convergence failure.

            if ( INFO != 0 ) {
               return;
            }
         } // 40
      } // 50

      }
