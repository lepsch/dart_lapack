      void zlaed7(final int N, final int CUTPNT, final int QSIZ, final int TLVLS, final int CURLVL, final int CURPBM, final int D, final Matrix<double> Q, final int LDQ, final int RHO, final int INDXQ, final int QSTORE, final int QPTR, final int PRMPTR, final int PERM, final int GIVPTR, final int GIVCOL, final int GIVNUM, final Array<double> _WORK, final Array<double> RWORK, final Array<int> IWORK, final Box<int> INFO,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                CURLVL, CURPBM, CUTPNT, INFO, LDQ, N, QSIZ, TLVLS;
      double             RHO;
      int                GIVCOL( 2, * ), GIVPTR( * ), INDXQ( * ), IWORK( * ), PERM( * ), PRMPTR( * ), QPTR( * );
      double             D( * ), GIVNUM( 2, * ), QSTORE( * ), RWORK( * );
      Complex         Q( LDQ, * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                COLTYP, CURR, I, IDLMDA, INDX, INDXC, INDXP, IQ, IW, IZ, K, N1, N2, PTR;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLAED9, DLAEDA, DLAMRG, XERBLA, ZLACRM, ZLAED8
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input parameters.

      INFO = 0;

      // IF( ICOMPQ < 0 || ICOMPQ > 1 ) THEN
      //    INFO = -1
      // ELSE IF( N < 0 ) THEN
      if ( N < 0 ) {
         INFO = -1;
      } else if ( min( 1, N ) > CUTPNT || N < CUTPNT ) {
         INFO = -2;
      } else if ( QSIZ < N ) {
         INFO = -3;
      } else if ( LDQ < max( 1, N ) ) {
         INFO = -9;
      }
      if ( INFO != 0 ) {
         xerbla('ZLAED7', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // The following values are for bookkeeping purposes only.  They are
      // integer pointers which indicate the portion of the workspace
      // used by a particular array in DLAED2 and SLAED3.

      IZ = 1;
      IDLMDA = IZ + N;
      IW = IDLMDA + N;
      IQ = IW + N;

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
      dlaeda(N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR, GIVCOL, GIVNUM, QSTORE, QPTR, RWORK( IZ ), RWORK( IZ+N ), INFO );

      // When solving the final problem, we no longer need the stored data,
      // so we will overwrite the data from this level onto the previously
      // used storage space.

      if ( CURLVL == TLVLS ) {
         QPTR[CURR] = 1;
         PRMPTR[CURR] = 1;
         GIVPTR[CURR] = 1;
      }

      // Sort and Deflate eigenvalues.

      zlaed8(K, N, QSIZ, Q, LDQ, D, RHO, CUTPNT, RWORK( IZ ), RWORK( IDLMDA ), WORK, QSIZ, RWORK( IW ), IWORK( INDXP ), IWORK( INDX ), INDXQ, PERM( PRMPTR( CURR ) ), GIVPTR( CURR+1 ), GIVCOL( 1, GIVPTR( CURR ) ), GIVNUM( 1, GIVPTR( CURR ) ), INFO );
      PRMPTR[CURR+1] = PRMPTR( CURR ) + N;
      GIVPTR[CURR+1] = GIVPTR( CURR+1 ) + GIVPTR( CURR );

      // Solve Secular Equation.

      if ( K != 0 ) {
         dlaed9(K, 1, K, N, D, RWORK( IQ ), K, RHO, RWORK( IDLMDA ), RWORK( IW ), QSTORE( QPTR( CURR ) ), K, INFO );
         zlacrm(QSIZ, K, WORK, QSIZ, QSTORE( QPTR( CURR ) ), K, Q, LDQ, RWORK( IQ ) );
         QPTR[CURR+1] = QPTR( CURR ) + K**2;
         if ( INFO != 0 ) {
            return;
         }

      // Prepare the INDXQ sorting permutation.

         N1 = K;
         N2 = N - K;
         dlamrg(N1, N2, D, 1, -1, INDXQ );
      } else {
         QPTR[CURR+1] = QPTR( CURR );
         for (I = 1; I <= N; I++) { // 20
            INDXQ[I] = I;
         } // 20
      }

      }
