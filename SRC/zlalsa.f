      void zlalsa(ICOMPQ, SMLSIZ, N, NRHS, B, LDB, BX, LDBX, U, LDU, VT, K, DIFL, DIFR, Z, POLES, GIVPTR, GIVCOL, LDGCOL, PERM, GIVNUM, C, S, RWORK, IWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                ICOMPQ, INFO, LDB, LDBX, LDGCOL, LDU, N, NRHS, SMLSIZ;
      // ..
      // .. Array Arguments ..
      int                GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ), K( * ), PERM( LDGCOL, * );
      double             C( * ), DIFL( LDU, * ), DIFR( LDU, * ), GIVNUM( LDU, * ), POLES( LDU, * ), RWORK( * ), S( * ), U( LDU, * ), VT( LDU, * ), Z( LDU, * );
      Complex         B( LDB, * ), BX( LDBX, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, I1, IC, IM1, INODE, J, JCOL, JIMAG, JREAL, JROW, LF, LL, LVL, LVL2, ND, NDB1, NDIML, NDIMR, NL, NLF, NLP1, NLVL, NR, NRF, NRP1, SQRE;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLASDT, XERBLA, ZCOPY, ZLALS0
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, DIMAG
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
         xerbla('ZLALSA', -INFO );
         return;
      }

      // Book-keeping and  setting up the computation tree.

      INODE = 1;
      NDIML = INODE + N;
      NDIMR = NDIML + N;

      dlasdt(N, NLVL, ND, IWORK( INODE ), IWORK( NDIML ), IWORK( NDIMR ), SMLSIZ );

      // The following code applies back the left singular vector factors.
      // For applying back the right singular vector factors, go to 170.

      if ( ICOMPQ == 1 ) {
         GO TO 170;
      }

      // The nodes on the bottom level of the tree were solved
      // by DLASDQ. The corresponding left and right singular vector
      // matrices are in explicit form. First apply back the left
      // singular vector matrices.

      NDB1 = ( ND+1 ) / 2;
      for (I = NDB1; I <= ND; I++) { // 130

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

         // Since B and BX are complex, the following call to DGEMM
         // is performed in two steps (real and imaginary parts).

         // CALL DGEMM( 'T', 'N', NL, NRHS, NL, ONE, U( NLF, 1 ), LDU,
      // $               B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX )

         J = NL*NRHS*2;
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 20
            for (JROW = NLF; JROW <= NLF + NL - 1; JROW++) { // 10
               J = J + 1;
               RWORK( J ) = DBLE( B( JROW, JCOL ) );
            } // 10
         } // 20
         dgemm('T', 'N', NL, NRHS, NL, ONE, U( NLF, 1 ), LDU, RWORK( 1+NL*NRHS*2 ), NL, ZERO, RWORK( 1 ), NL );
         J = NL*NRHS*2;
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 40
            for (JROW = NLF; JROW <= NLF + NL - 1; JROW++) { // 30
               J = J + 1;
               RWORK( J ) = DIMAG( B( JROW, JCOL ) );
            } // 30
         } // 40
         dgemm('T', 'N', NL, NRHS, NL, ONE, U( NLF, 1 ), LDU, RWORK( 1+NL*NRHS*2 ), NL, ZERO, RWORK( 1+NL*NRHS ), NL );
         JREAL = 0;
         JIMAG = NL*NRHS;
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 60
            for (JROW = NLF; JROW <= NLF + NL - 1; JROW++) { // 50
               JREAL = JREAL + 1;
               JIMAG = JIMAG + 1;
               BX( JROW, JCOL ) = DCMPLX( RWORK( JREAL ), RWORK( JIMAG ) );
            } // 50
         } // 60

         // Since B and BX are complex, the following call to DGEMM
         // is performed in two steps (real and imaginary parts).

         // CALL DGEMM( 'T', 'N', NR, NRHS, NR, ONE, U( NRF, 1 ), LDU,
// $               B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX )

         J = NR*NRHS*2;
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 80
            for (JROW = NRF; JROW <= NRF + NR - 1; JROW++) { // 70
               J = J + 1;
               RWORK( J ) = DBLE( B( JROW, JCOL ) );
            } // 70
         } // 80
         dgemm('T', 'N', NR, NRHS, NR, ONE, U( NRF, 1 ), LDU, RWORK( 1+NR*NRHS*2 ), NR, ZERO, RWORK( 1 ), NR );
         J = NR*NRHS*2;
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 100
            for (JROW = NRF; JROW <= NRF + NR - 1; JROW++) { // 90
               J = J + 1;
               RWORK( J ) = DIMAG( B( JROW, JCOL ) );
            } // 90
         } // 100
         dgemm('T', 'N', NR, NRHS, NR, ONE, U( NRF, 1 ), LDU, RWORK( 1+NR*NRHS*2 ), NR, ZERO, RWORK( 1+NR*NRHS ), NR );
         JREAL = 0;
         JIMAG = NR*NRHS;
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 120
            for (JROW = NRF; JROW <= NRF + NR - 1; JROW++) { // 110
               JREAL = JREAL + 1;
               JIMAG = JIMAG + 1;
               BX( JROW, JCOL ) = DCMPLX( RWORK( JREAL ), RWORK( JIMAG ) );
            } // 110
         } // 120

      } // 130

      // Next copy the rows of B that correspond to unchanged rows
      // in the bidiagonal matrix to BX.

      for (I = 1; I <= ND; I++) { // 140
         IC = IWORK( INODE+I-1 );
         zcopy(NRHS, B( IC, 1 ), LDB, BX( IC, 1 ), LDBX );
      } // 140

      // Finally go through the left singular vector matrices of all
      // the other subproblems bottom-up on the tree.

      J = 2**NLVL;
      SQRE = 0;

      DO 160 LVL = NLVL, 1, -1;
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
         for (I = LF; I <= LL; I++) { // 150
            IM1 = I - 1;
            IC = IWORK( INODE+IM1 );
            NL = IWORK( NDIML+IM1 );
            NR = IWORK( NDIMR+IM1 );
            NLF = IC - NL;
            NRF = IC + 1;
            J = J - 1;
            zlals0(ICOMPQ, NL, NR, SQRE, NRHS, BX( NLF, 1 ), LDBX, B( NLF, 1 ), LDB, PERM( NLF, LVL ), GIVPTR( J ), GIVCOL( NLF, LVL2 ), LDGCOL, GIVNUM( NLF, LVL2 ), LDU, POLES( NLF, LVL2 ), DIFL( NLF, LVL ), DIFR( NLF, LVL2 ), Z( NLF, LVL ), K( J ), C( J ), S( J ), RWORK, INFO );
         } // 150
      } // 160
      GO TO 330;

      // ICOMPQ = 1: applying back the right singular vector factors.

      } // 170

      // First now go through the right singular vector matrices of all
      // the tree nodes top-down.

      J = 0;
      for (LVL = 1; LVL <= NLVL; LVL++) { // 190
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
         DO 180 I = LL, LF, -1;
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
            zlals0(ICOMPQ, NL, NR, SQRE, NRHS, B( NLF, 1 ), LDB, BX( NLF, 1 ), LDBX, PERM( NLF, LVL ), GIVPTR( J ), GIVCOL( NLF, LVL2 ), LDGCOL, GIVNUM( NLF, LVL2 ), LDU, POLES( NLF, LVL2 ), DIFL( NLF, LVL ), DIFR( NLF, LVL2 ), Z( NLF, LVL ), K( J ), C( J ), S( J ), RWORK, INFO );
         } // 180
      } // 190

      // The nodes on the bottom level of the tree were solved
      // by DLASDQ. The corresponding right singular vector
      // matrices are in explicit form. Apply them back.

      NDB1 = ( ND+1 ) / 2;
      for (I = NDB1; I <= ND; I++) { // 320
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

         // Since B and BX are complex, the following call to DGEMM is
         // performed in two steps (real and imaginary parts).

         // CALL DGEMM( 'T', 'N', NLP1, NRHS, NLP1, ONE, VT( NLF, 1 ), LDU,
// $               B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX )

         J = NLP1*NRHS*2;
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 210
            for (JROW = NLF; JROW <= NLF + NLP1 - 1; JROW++) { // 200
               J = J + 1;
               RWORK( J ) = DBLE( B( JROW, JCOL ) );
            } // 200
         } // 210
         dgemm('T', 'N', NLP1, NRHS, NLP1, ONE, VT( NLF, 1 ), LDU, RWORK( 1+NLP1*NRHS*2 ), NLP1, ZERO, RWORK( 1 ), NLP1 );
         J = NLP1*NRHS*2;
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 230
            for (JROW = NLF; JROW <= NLF + NLP1 - 1; JROW++) { // 220
               J = J + 1;
               RWORK( J ) = DIMAG( B( JROW, JCOL ) );
            } // 220
         } // 230
         dgemm('T', 'N', NLP1, NRHS, NLP1, ONE, VT( NLF, 1 ), LDU, RWORK( 1+NLP1*NRHS*2 ), NLP1, ZERO, RWORK( 1+NLP1*NRHS ), NLP1 );
         JREAL = 0;
         JIMAG = NLP1*NRHS;
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 250
            for (JROW = NLF; JROW <= NLF + NLP1 - 1; JROW++) { // 240
               JREAL = JREAL + 1;
               JIMAG = JIMAG + 1;
               BX( JROW, JCOL ) = DCMPLX( RWORK( JREAL ), RWORK( JIMAG ) );
            } // 240
         } // 250

         // Since B and BX are complex, the following call to DGEMM is
         // performed in two steps (real and imaginary parts).

         // CALL DGEMM( 'T', 'N', NRP1, NRHS, NRP1, ONE, VT( NRF, 1 ), LDU,
// $               B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX )

         J = NRP1*NRHS*2;
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 270
            for (JROW = NRF; JROW <= NRF + NRP1 - 1; JROW++) { // 260
               J = J + 1;
               RWORK( J ) = DBLE( B( JROW, JCOL ) );
            } // 260
         } // 270
         dgemm('T', 'N', NRP1, NRHS, NRP1, ONE, VT( NRF, 1 ), LDU, RWORK( 1+NRP1*NRHS*2 ), NRP1, ZERO, RWORK( 1 ), NRP1 );
         J = NRP1*NRHS*2;
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 290
            for (JROW = NRF; JROW <= NRF + NRP1 - 1; JROW++) { // 280
               J = J + 1;
               RWORK( J ) = DIMAG( B( JROW, JCOL ) );
            } // 280
         } // 290
         dgemm('T', 'N', NRP1, NRHS, NRP1, ONE, VT( NRF, 1 ), LDU, RWORK( 1+NRP1*NRHS*2 ), NRP1, ZERO, RWORK( 1+NRP1*NRHS ), NRP1 );
         JREAL = 0;
         JIMAG = NRP1*NRHS;
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 310
            for (JROW = NRF; JROW <= NRF + NRP1 - 1; JROW++) { // 300
               JREAL = JREAL + 1;
               JIMAG = JIMAG + 1;
               BX( JROW, JCOL ) = DCMPLX( RWORK( JREAL ), RWORK( JIMAG ) );
            } // 300
         } // 310

      } // 320

      } // 330

      return;

      // End of ZLALSA

      }
