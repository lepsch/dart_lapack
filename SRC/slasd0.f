      SUBROUTINE SLASD0( N, SQRE, D, E, U, LDU, VT, LDVT, SMLSIZ, IWORK, WORK, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDU, LDVT, N, SMLSIZ, SQRE;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               D( * ), E( * ), U( LDU, * ), VT( LDVT, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, I1, IC, IDXQ, IDXQC, IM1, INODE, ITEMP, IWK, J, LF, LL, LVL, M, NCC, ND, NDB1, NDIML, NDIMR, NL, NLF, NLP1, NLVL, NR, NRF, NRP1, SQREI;
      REAL               ALPHA, BETA
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASD1, SLASDQ, SLASDT, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0

      if ( N.LT.0 ) {
         INFO = -1
      } else if ( ( SQRE.LT.0 ) .OR. ( SQRE.GT.1 ) ) {
         INFO = -2
      }

      M = N + SQRE

      if ( LDU.LT.N ) {
         INFO = -6
      } else if ( LDVT.LT.M ) {
         INFO = -8
      } else if ( SMLSIZ.LT.3 ) {
         INFO = -9
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'SLASD0', -INFO )
         RETURN
      }

      // If the input matrix is too small, call SLASDQ to find the SVD.

      if ( N.LE.SMLSIZ ) {
         CALL SLASDQ( 'U', SQRE, N, M, N, 0, D, E, VT, LDVT, U, LDU, U, LDU, WORK, INFO )
         RETURN
      }

      // Set up the computation tree.

      INODE = 1
      NDIML = INODE + N
      NDIMR = NDIML + N
      IDXQ = NDIMR + N
      IWK = IDXQ + N
      CALL SLASDT( N, NLVL, ND, IWORK( INODE ), IWORK( NDIML ), IWORK( NDIMR ), SMLSIZ )

      // For the nodes on bottom level of the tree, solve
     t // heir subproblems by SLASDQ.

      NDB1 = ( ND+1 ) / 2
      NCC = 0
      DO 30 I = NDB1, ND

      // IC : center row of each node
      // NL : number of rows of left  subproblem
      // NR : number of rows of right subproblem
      // NLF: starting row of the left   subproblem
      // NRF: starting row of the right  subproblem

         I1 = I - 1
         IC = IWORK( INODE+I1 )
         NL = IWORK( NDIML+I1 )
         NLP1 = NL + 1
         NR = IWORK( NDIMR+I1 )
         NRP1 = NR + 1
         NLF = IC - NL
         NRF = IC + 1
         SQREI = 1
         CALL SLASDQ( 'U', SQREI, NL, NLP1, NL, NCC, D( NLF ), E( NLF ), VT( NLF, NLF ), LDVT, U( NLF, NLF ), LDU, U( NLF, NLF ), LDU, WORK, INFO )
         if ( INFO.NE.0 ) {
            RETURN
         }
         ITEMP = IDXQ + NLF - 2
         DO 10 J = 1, NL
            IWORK( ITEMP+J ) = J
   10    CONTINUE
         if ( I.EQ.ND ) {
            SQREI = SQRE
         } else {
            SQREI = 1
         }
         NRP1 = NR + SQREI
         CALL SLASDQ( 'U', SQREI, NR, NRP1, NR, NCC, D( NRF ), E( NRF ), VT( NRF, NRF ), LDVT, U( NRF, NRF ), LDU, U( NRF, NRF ), LDU, WORK, INFO )
         if ( INFO.NE.0 ) {
            RETURN
         }
         ITEMP = IDXQ + IC
         DO 20 J = 1, NR
            IWORK( ITEMP+J-1 ) = J
   20    CONTINUE
   30 CONTINUE

      // Now conquer each subproblem bottom-up.

      DO 50 LVL = NLVL, 1, -1

         // Find the first node LF and last node LL on the
         // current level LVL.

         if ( LVL.EQ.1 ) {
            LF = 1
            LL = 1
         } else {
            LF = 2**( LVL-1 )
            LL = 2*LF - 1
         }
         DO 40 I = LF, LL
            IM1 = I - 1
            IC = IWORK( INODE+IM1 )
            NL = IWORK( NDIML+IM1 )
            NR = IWORK( NDIMR+IM1 )
            NLF = IC - NL
            if ( ( SQRE.EQ.0 ) .AND. ( I.EQ.LL ) ) {
               SQREI = SQRE
            } else {
               SQREI = 1
            }
            IDXQC = IDXQ + NLF - 1
            ALPHA = D( IC )
            BETA = E( IC )
            CALL SLASD1( NL, NR, SQREI, D( NLF ), ALPHA, BETA, U( NLF, NLF ), LDU, VT( NLF, NLF ), LDVT, IWORK( IDXQC ), IWORK( IWK ), WORK, INFO )

      // Report the possible convergence failure.

            if ( INFO.NE.0 ) {
               RETURN
            }
   40    CONTINUE
   50 CONTINUE

      RETURN

      // End of SLASD0

      }
