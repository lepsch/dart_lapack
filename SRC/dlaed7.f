      void dlaed7(ICOMPQ, N, QSIZ, TLVLS, CURLVL, CURPBM, D, Q, LDQ, INDXQ, RHO, CUTPNT, QSTORE, QPTR, PRMPTR, PERM, GIVPTR, GIVCOL, GIVNUM, WORK, IWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                CURLVL, CURPBM, CUTPNT, ICOMPQ, INFO, LDQ, N, QSIZ, TLVLS;
      double             RHO;
      // ..
      // .. Array Arguments ..
      int                GIVCOL( 2, * ), GIVPTR( * ), INDXQ( * ), IWORK( * ), PERM( * ), PRMPTR( * ), QPTR( * );
      double             D( * ), GIVNUM( 2, * ), Q( LDQ, * ), QSTORE( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                COLTYP, CURR, I, IDLMDA, INDX, INDXC, INDXP, IQ2, IS, IW, IZ, K, LDQ2, N1, N2, PTR;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLAED8, DLAED9, DLAEDA, DLAMRG, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;

      if ( ICOMPQ < 0 || ICOMPQ > 1 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( ICOMPQ == 1 && QSIZ < N ) {
         INFO = -3;
      } else if ( LDQ < max( 1, N ) ) {
         INFO = -9;
      } else if ( min( 1, N ) > CUTPNT || N < CUTPNT ) {
         INFO = -12;
      }
      if ( INFO != 0 ) {
         xerbla('DLAED7', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // The following values are for bookkeeping purposes only.  They are
      // integer pointers which indicate the portion of the workspace
      // used by a particular array in DLAED8 and DLAED9.

      if ( ICOMPQ == 1 ) {
         LDQ2 = QSIZ;
      } else {
         LDQ2 = N;
      }

      IZ = 1;
      IDLMDA = IZ + N;
      IW = IDLMDA + N;
      IQ2 = IW + N;
      IS = IQ2 + N*LDQ2;

      INDX = 1;
      INDXC = INDX + N;
      COLTYP = INDXC + N;
      INDXP = COLTYP + N;

      // Form the z-vector which consists of the last row of Q_1 and the
      // first row of Q_2.

      PTR = 1 + 2**TLVLS;
      for (I = 1; I <= CURLVL - 1; I++) { // 10
         PTR = PTR + 2**( TLVLS-I );
      } // 10
      CURR = PTR + CURPBM;
      dlaeda(N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR, GIVCOL, GIVNUM, QSTORE, QPTR, WORK( IZ ), WORK( IZ+N ), INFO );

      // When solving the final problem, we no longer need the stored data,
      // so we will overwrite the data from this level onto the previously
      // used storage space.

      if ( CURLVL == TLVLS ) {
         QPTR( CURR ) = 1;
         PRMPTR( CURR ) = 1;
         GIVPTR( CURR ) = 1;
      }

      // Sort and Deflate eigenvalues.

      dlaed8(ICOMPQ, K, N, QSIZ, D, Q, LDQ, INDXQ, RHO, CUTPNT, WORK( IZ ), WORK( IDLMDA ), WORK( IQ2 ), LDQ2, WORK( IW ), PERM( PRMPTR( CURR ) ), GIVPTR( CURR+1 ), GIVCOL( 1, GIVPTR( CURR ) ), GIVNUM( 1, GIVPTR( CURR ) ), IWORK( INDXP ), IWORK( INDX ), INFO );
      PRMPTR( CURR+1 ) = PRMPTR( CURR ) + N;
      GIVPTR( CURR+1 ) = GIVPTR( CURR+1 ) + GIVPTR( CURR );

      // Solve Secular Equation.

      if ( K != 0 ) {
         CALL DLAED9( K, 1, K, N, D, WORK( IS ), K, RHO, WORK( IDLMDA ), WORK( IW ), QSTORE( QPTR( CURR ) ), K, INFO )          IF( INFO != 0 ) GO TO 30;
         if ( ICOMPQ == 1 ) {
            dgemm('N', 'N', QSIZ, K, K, ONE, WORK( IQ2 ), LDQ2, QSTORE( QPTR( CURR ) ), K, ZERO, Q, LDQ );
         }
         QPTR( CURR+1 ) = QPTR( CURR ) + K**2;

      // Prepare the INDXQ sorting permutation.

         N1 = K;
         N2 = N - K;
         dlamrg(N1, N2, D, 1, -1, INDXQ );
      } else {
         QPTR( CURR+1 ) = QPTR( CURR );
         for (I = 1; I <= N; I++) { // 20
            INDXQ( I ) = I;
         } // 20
      }

      } // 30
      return;
      }
