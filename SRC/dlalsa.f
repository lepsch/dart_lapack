      SUBROUTINE DLALSA( ICOMPQ, SMLSIZ, N, NRHS, B, LDB, BX, LDBX, U, LDU, VT, K, DIFL, DIFR, Z, POLES, GIVPTR, GIVCOL, LDGCOL, PERM, GIVNUM, C, S, WORK, IWORK, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                ICOMPQ, INFO, LDB, LDBX, LDGCOL, LDU, N, NRHS, SMLSIZ;
      // ..
      // .. Array Arguments ..
      int                GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ), K( * ), PERM( LDGCOL, * );
      double             B( LDB, * ), BX( LDBX, * ), C( * ), DIFL( LDU, * ), DIFR( LDU, * ), GIVNUM( LDU, * ), POLES( LDU, * ), S( * ), U( LDU, * ), VT( LDU, * ), WORK( * ), Z( LDU, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, I1, IC, IM1, INODE, J, LF, LL, LVL, LVL2, ND, NDB1, NDIML, NDIMR, NL, NLF, NLP1, NLVL, NR, NRF, NRP1, SQRE;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGEMM, DLALS0, DLASDT, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;

      if ( ( ICOMPQ < 0 ) || ( ICOMPQ > 1 ) ) {
         INFO = -1;
      } else if ( SMLSIZ < 3 ) {
         INFO = -2;
      } else if ( N < SMLSIZ ) {
         INFO = -3;
      } else if ( NRHS < 1 ) {
         INFO = -4;
      } else if ( LDB < N ) {
         INFO = -6;
      } else if ( LDBX < N ) {
         INFO = -8;
      } else if ( LDU < N ) {
         INFO = -10;
      } else if ( LDGCOL < N ) {
         INFO = -19;
      }
      if ( INFO != 0 ) {
         xerbla('DLALSA', -INFO );
         return;
      }

      // Book-keeping and  setting up the computation tree.

      INODE = 1;
      NDIML = INODE + N;
      NDIMR = NDIML + N;

      dlasdt(N, NLVL, ND, IWORK( INODE ), IWORK( NDIML ), IWORK( NDIMR ), SMLSIZ );

      // The following code applies back the left singular vector factors.
      // For applying back the right singular vector factors, go to 50.

      if ( ICOMPQ == 1 ) {
         GO TO 50;
      }

      // The nodes on the bottom level of the tree were solved
      // by DLASDQ. The corresponding left and right singular vector
      // matrices are in explicit form. First apply back the left
      // singular vector matrices.

      NDB1 = ( ND+1 ) / 2;
      for (I = NDB1; I <= ND; I++) { // 10

         // IC : center row of each node
         // NL : number of rows of left  subproblem
         // NR : number of rows of right subproblem
         // NLF: starting row of the left   subproblem
         // NRF: starting row of the right  subproblem

         I1 = I - 1;
         IC = IWORK( INODE+I1 );
         NL = IWORK( NDIML+I1 );
         NR = IWORK( NDIMR+I1 );
         NLF = IC - NL;
         NRF = IC + 1;
         dgemm('T', 'N', NL, NRHS, NL, ONE, U( NLF, 1 ), LDU, B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX );
         dgemm('T', 'N', NR, NRHS, NR, ONE, U( NRF, 1 ), LDU, B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX );
      } // 10

      // Next copy the rows of B that correspond to unchanged rows
      // in the bidiagonal matrix to BX.

      for (I = 1; I <= ND; I++) { // 20
         IC = IWORK( INODE+I-1 );
         dcopy(NRHS, B( IC, 1 ), LDB, BX( IC, 1 ), LDBX );
      } // 20

      // Finally go through the left singular vector matrices of all
      // the other subproblems bottom-up on the tree.

      J = 2**NLVL;
      SQRE = 0;

      DO 40 LVL = NLVL, 1, -1;
         LVL2 = 2*LVL - 1;

         // find the first node LF and last node LL on
         // the current level LVL

         if ( LVL == 1 ) {
            LF = 1;
            LL = 1;
         } else {
            LF = 2**( LVL-1 );
            LL = 2*LF - 1;
         }
         for (I = LF; I <= LL; I++) { // 30
            IM1 = I - 1;
            IC = IWORK( INODE+IM1 );
            NL = IWORK( NDIML+IM1 );
            NR = IWORK( NDIMR+IM1 );
            NLF = IC - NL;
            NRF = IC + 1;
            J = J - 1;
            dlals0(ICOMPQ, NL, NR, SQRE, NRHS, BX( NLF, 1 ), LDBX, B( NLF, 1 ), LDB, PERM( NLF, LVL ), GIVPTR( J ), GIVCOL( NLF, LVL2 ), LDGCOL, GIVNUM( NLF, LVL2 ), LDU, POLES( NLF, LVL2 ), DIFL( NLF, LVL ), DIFR( NLF, LVL2 ), Z( NLF, LVL ), K( J ), C( J ), S( J ), WORK, INFO );
         } // 30
      } // 40
      GO TO 90;

      // ICOMPQ = 1: applying back the right singular vector factors.

      } // 50

      // First now go through the right singular vector matrices of all
      // the tree nodes top-down.

      J = 0;
      for (LVL = 1; LVL <= NLVL; LVL++) { // 70
         LVL2 = 2*LVL - 1;

         // Find the first node LF and last node LL on
         // the current level LVL.

         if ( LVL == 1 ) {
            LF = 1;
            LL = 1;
         } else {
            LF = 2**( LVL-1 );
            LL = 2*LF - 1;
         }
         DO 60 I = LL, LF, -1;
            IM1 = I - 1;
            IC = IWORK( INODE+IM1 );
            NL = IWORK( NDIML+IM1 );
            NR = IWORK( NDIMR+IM1 );
            NLF = IC - NL;
            NRF = IC + 1;
            if ( I == LL ) {
               SQRE = 0;
            } else {
               SQRE = 1;
            }
            J = J + 1;
            dlals0(ICOMPQ, NL, NR, SQRE, NRHS, B( NLF, 1 ), LDB, BX( NLF, 1 ), LDBX, PERM( NLF, LVL ), GIVPTR( J ), GIVCOL( NLF, LVL2 ), LDGCOL, GIVNUM( NLF, LVL2 ), LDU, POLES( NLF, LVL2 ), DIFL( NLF, LVL ), DIFR( NLF, LVL2 ), Z( NLF, LVL ), K( J ), C( J ), S( J ), WORK, INFO );
         } // 60
      } // 70

      // The nodes on the bottom level of the tree were solved
      // by DLASDQ. The corresponding right singular vector
      // matrices are in explicit form. Apply them back.

      NDB1 = ( ND+1 ) / 2;
      for (I = NDB1; I <= ND; I++) { // 80
         I1 = I - 1;
         IC = IWORK( INODE+I1 );
         NL = IWORK( NDIML+I1 );
         NR = IWORK( NDIMR+I1 );
         NLP1 = NL + 1;
         if ( I == ND ) {
            NRP1 = NR;
         } else {
            NRP1 = NR + 1;
         }
         NLF = IC - NL;
         NRF = IC + 1;
         dgemm('T', 'N', NLP1, NRHS, NLP1, ONE, VT( NLF, 1 ), LDU, B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX );
         dgemm('T', 'N', NRP1, NRHS, NRP1, ONE, VT( NRF, 1 ), LDU, B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX );
      } // 80

      } // 90

      return;

      // End of DLALSA

      }
