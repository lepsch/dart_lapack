      SUBROUTINE CLALSA( ICOMPQ, SMLSIZ, N, NRHS, B, LDB, BX, LDBX, U, LDU, VT, K, DIFL, DIFR, Z, POLES, GIVPTR, GIVCOL, LDGCOL, PERM, GIVNUM, C, S, RWORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                ICOMPQ, INFO, LDB, LDBX, LDGCOL, LDU, N, NRHS, SMLSIZ;
      // ..
      // .. Array Arguments ..
      int                GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ), K( * ), PERM( LDGCOL, * );
      REAL               C( * ), DIFL( LDU, * ), DIFR( LDU, * ), GIVNUM( LDU, * ), POLES( LDU, * ), RWORK( * ), S( * ), U( LDU, * ), VT( LDU, * ), Z( LDU, * );
      COMPLEX            B( LDB, * ), BX( LDBX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      // ..
      // .. Local Scalars ..
      int                I, I1, IC, IM1, INODE, J, JCOL, JIMAG, JREAL, JROW, LF, LL, LVL, LVL2, ND, NDB1, NDIML, NDIMR, NL, NLF, NLP1, NLVL, NR, NRF, NRP1, SQRE;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CLALS0, SGEMM, SLASDT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, CMPLX, REAL
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
         xerbla('CLALSA', -INFO );
         RETURN
      }

      // Book-keeping and  setting up the computation tree.

      INODE = 1
      NDIML = INODE + N
      NDIMR = NDIML + N

      slasdt(N, NLVL, ND, IWORK( INODE ), IWORK( NDIML ), IWORK( NDIMR ), SMLSIZ );

      // The following code applies back the left singular vector factors.
      // For applying back the right singular vector factors, go to 170.

      if ( ICOMPQ.EQ.1 ) {
         GO TO 170
      }

      // The nodes on the bottom level of the tree were solved
      // by SLASDQ. The corresponding left and right singular vector
      // matrices are in explicit form. First apply back the left
      // singular vector matrices.

      NDB1 = ( ND+1 ) / 2
      for (I = NDB1; I <= ND; I++) { // 130

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

         // Since B and BX are complex, the following call to SGEMM
         // is performed in two steps (real and imaginary parts).

         // CALL SGEMM( 'T', 'N', NL, NRHS, NL, ONE, U( NLF, 1 ), LDU,
      // $               B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX )

         J = NL*NRHS*2
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 20
            DO 10 JROW = NLF, NLF + NL - 1
               J = J + 1
               RWORK( J ) = REAL( B( JROW, JCOL ) )
   10       CONTINUE
   20    CONTINUE
         sgemm('T', 'N', NL, NRHS, NL, ONE, U( NLF, 1 ), LDU, RWORK( 1+NL*NRHS*2 ), NL, ZERO, RWORK( 1 ), NL );
         J = NL*NRHS*2
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 40
            DO 30 JROW = NLF, NLF + NL - 1
               J = J + 1
               RWORK( J ) = AIMAG( B( JROW, JCOL ) )
   30       CONTINUE
   40    CONTINUE
         sgemm('T', 'N', NL, NRHS, NL, ONE, U( NLF, 1 ), LDU, RWORK( 1+NL*NRHS*2 ), NL, ZERO, RWORK( 1+NL*NRHS ), NL );
         JREAL = 0
         JIMAG = NL*NRHS
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 60
            DO 50 JROW = NLF, NLF + NL - 1
               JREAL = JREAL + 1
               JIMAG = JIMAG + 1
               BX( JROW, JCOL ) = CMPLX( RWORK( JREAL ), RWORK( JIMAG ) )
   50       CONTINUE
   60    CONTINUE

         // Since B and BX are complex, the following call to SGEMM
         // is performed in two steps (real and imaginary parts).

         // CALL SGEMM( 'T', 'N', NR, NRHS, NR, ONE, U( NRF, 1 ), LDU,
*    $               B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX )

         J = NR*NRHS*2
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 80
            DO 70 JROW = NRF, NRF + NR - 1
               J = J + 1
               RWORK( J ) = REAL( B( JROW, JCOL ) )
   70       CONTINUE
   80    CONTINUE
         sgemm('T', 'N', NR, NRHS, NR, ONE, U( NRF, 1 ), LDU, RWORK( 1+NR*NRHS*2 ), NR, ZERO, RWORK( 1 ), NR );
         J = NR*NRHS*2
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 100
            DO 90 JROW = NRF, NRF + NR - 1
               J = J + 1
               RWORK( J ) = AIMAG( B( JROW, JCOL ) )
   90       CONTINUE
  100    CONTINUE
         sgemm('T', 'N', NR, NRHS, NR, ONE, U( NRF, 1 ), LDU, RWORK( 1+NR*NRHS*2 ), NR, ZERO, RWORK( 1+NR*NRHS ), NR );
         JREAL = 0
         JIMAG = NR*NRHS
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 120
            DO 110 JROW = NRF, NRF + NR - 1
               JREAL = JREAL + 1
               JIMAG = JIMAG + 1
               BX( JROW, JCOL ) = CMPLX( RWORK( JREAL ), RWORK( JIMAG ) )
  110       CONTINUE
  120    CONTINUE

  130 CONTINUE

      // Next copy the rows of B that correspond to unchanged rows
      // in the bidiagonal matrix to BX.

      for (I = 1; I <= ND; I++) { // 140
         IC = IWORK( INODE+I-1 )
         ccopy(NRHS, B( IC, 1 ), LDB, BX( IC, 1 ), LDBX );
  140 CONTINUE

      // Finally go through the left singular vector matrices of all
      // the other subproblems bottom-up on the tree.

      J = 2**NLVL
      SQRE = 0

      DO 160 LVL = NLVL, 1, -1
         LVL2 = 2*LVL - 1

         // find the first node LF and last node LL on
         // the current level LVL

         if ( LVL.EQ.1 ) {
            LF = 1
            LL = 1
         } else {
            LF = 2**( LVL-1 )
            LL = 2*LF - 1
         }
         for (I = LF; I <= LL; I++) { // 150
            IM1 = I - 1
            IC = IWORK( INODE+IM1 )
            NL = IWORK( NDIML+IM1 )
            NR = IWORK( NDIMR+IM1 )
            NLF = IC - NL
            NRF = IC + 1
            J = J - 1
            clals0(ICOMPQ, NL, NR, SQRE, NRHS, BX( NLF, 1 ), LDBX, B( NLF, 1 ), LDB, PERM( NLF, LVL ), GIVPTR( J ), GIVCOL( NLF, LVL2 ), LDGCOL, GIVNUM( NLF, LVL2 ), LDU, POLES( NLF, LVL2 ), DIFL( NLF, LVL ), DIFR( NLF, LVL2 ), Z( NLF, LVL ), K( J ), C( J ), S( J ), RWORK, INFO );
  150    CONTINUE
  160 CONTINUE
      GO TO 330

      // ICOMPQ = 1: applying back the right singular vector factors.

  170 CONTINUE

      // First now go through the right singular vector matrices of all
      // the tree nodes top-down.

      J = 0
      for (LVL = 1; LVL <= NLVL; LVL++) { // 190
         LVL2 = 2*LVL - 1

         // Find the first node LF and last node LL on
         // the current level LVL.

         if ( LVL.EQ.1 ) {
            LF = 1
            LL = 1
         } else {
            LF = 2**( LVL-1 )
            LL = 2*LF - 1
         }
         DO 180 I = LL, LF, -1
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
            clals0(ICOMPQ, NL, NR, SQRE, NRHS, B( NLF, 1 ), LDB, BX( NLF, 1 ), LDBX, PERM( NLF, LVL ), GIVPTR( J ), GIVCOL( NLF, LVL2 ), LDGCOL, GIVNUM( NLF, LVL2 ), LDU, POLES( NLF, LVL2 ), DIFL( NLF, LVL ), DIFR( NLF, LVL2 ), Z( NLF, LVL ), K( J ), C( J ), S( J ), RWORK, INFO );
  180    CONTINUE
  190 CONTINUE

      // The nodes on the bottom level of the tree were solved
      // by SLASDQ. The corresponding right singular vector
      // matrices are in explicit form. Apply them back.

      NDB1 = ( ND+1 ) / 2
      for (I = NDB1; I <= ND; I++) { // 320
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

         // Since B and BX are complex, the following call to SGEMM is
         // performed in two steps (real and imaginary parts).

         // CALL SGEMM( 'T', 'N', NLP1, NRHS, NLP1, ONE, VT( NLF, 1 ), LDU,
*    $               B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX )

         J = NLP1*NRHS*2
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 210
            DO 200 JROW = NLF, NLF + NLP1 - 1
               J = J + 1
               RWORK( J ) = REAL( B( JROW, JCOL ) )
  200       CONTINUE
  210    CONTINUE
         sgemm('T', 'N', NLP1, NRHS, NLP1, ONE, VT( NLF, 1 ), LDU, RWORK( 1+NLP1*NRHS*2 ), NLP1, ZERO, RWORK( 1 ), NLP1 );
         J = NLP1*NRHS*2
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 230
            DO 220 JROW = NLF, NLF + NLP1 - 1
               J = J + 1
               RWORK( J ) = AIMAG( B( JROW, JCOL ) )
  220       CONTINUE
  230    CONTINUE
         sgemm('T', 'N', NLP1, NRHS, NLP1, ONE, VT( NLF, 1 ), LDU, RWORK( 1+NLP1*NRHS*2 ), NLP1, ZERO, RWORK( 1+NLP1*NRHS ), NLP1 );
         JREAL = 0
         JIMAG = NLP1*NRHS
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 250
            DO 240 JROW = NLF, NLF + NLP1 - 1
               JREAL = JREAL + 1
               JIMAG = JIMAG + 1
               BX( JROW, JCOL ) = CMPLX( RWORK( JREAL ), RWORK( JIMAG ) )
  240       CONTINUE
  250    CONTINUE

         // Since B and BX are complex, the following call to SGEMM is
         // performed in two steps (real and imaginary parts).

         // CALL SGEMM( 'T', 'N', NRP1, NRHS, NRP1, ONE, VT( NRF, 1 ), LDU,
*    $               B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX )

         J = NRP1*NRHS*2
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 270
            DO 260 JROW = NRF, NRF + NRP1 - 1
               J = J + 1
               RWORK( J ) = REAL( B( JROW, JCOL ) )
  260       CONTINUE
  270    CONTINUE
         sgemm('T', 'N', NRP1, NRHS, NRP1, ONE, VT( NRF, 1 ), LDU, RWORK( 1+NRP1*NRHS*2 ), NRP1, ZERO, RWORK( 1 ), NRP1 );
         J = NRP1*NRHS*2
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 290
            DO 280 JROW = NRF, NRF + NRP1 - 1
               J = J + 1
               RWORK( J ) = AIMAG( B( JROW, JCOL ) )
  280       CONTINUE
  290    CONTINUE
         sgemm('T', 'N', NRP1, NRHS, NRP1, ONE, VT( NRF, 1 ), LDU, RWORK( 1+NRP1*NRHS*2 ), NRP1, ZERO, RWORK( 1+NRP1*NRHS ), NRP1 );
         JREAL = 0
         JIMAG = NRP1*NRHS
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 310
            DO 300 JROW = NRF, NRF + NRP1 - 1
               JREAL = JREAL + 1
               JIMAG = JIMAG + 1
               BX( JROW, JCOL ) = CMPLX( RWORK( JREAL ), RWORK( JIMAG ) )
  300       CONTINUE
  310    CONTINUE

  320 CONTINUE

  330 CONTINUE

      RETURN

      // End of CLALSA

      }
