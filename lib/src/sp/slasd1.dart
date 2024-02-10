      void slasd1(NL, NR, SQRE, D, ALPHA, BETA, final Matrix<double> U, final int LDU, final Matrix<double> VT, final int LDVT, IDXQ, final Array<int> IWORK, final Array<double> _WORK, final Box<int> INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDU, LDVT, NL, NR, SQRE;
      double               ALPHA, BETA;
      int                IDXQ( * ), IWORK( * );
      double               D( * ), U( LDU, * ), VT( LDVT, * ), WORK( * );
      // ..


      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                COLTYP, I, IDX, IDXC, IDXP, IQ, ISIGMA, IU2, IVT2, IZ, K, LDQ, LDU2, LDVT2, M, N, N1, N2;
      double               ORGNRM;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLAMRG, SLASCL, SLASD2, SLASD3, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX

      // Test the input parameters.

      INFO = 0;

      if ( NL < 1 ) {
         INFO = -1;
      } else if ( NR < 1 ) {
         INFO = -2;
      } else if ( ( SQRE < 0 ) || ( SQRE > 1 ) ) {
         INFO = -3;
      }
      if ( INFO != 0 ) {
         xerbla('SLASD1', -INFO );
         return;
      }

      N = NL + NR + 1;
      M = N + SQRE;

      // The following values are for bookkeeping purposes only.  They are
      // integer pointers which indicate the portion of the workspace
      // used by a particular array in SLASD2 and SLASD3.

      LDU2 = N;
      LDVT2 = M;

      IZ = 1;
      ISIGMA = IZ + M;
      IU2 = ISIGMA + N;
      IVT2 = IU2 + LDU2*N;
      IQ = IVT2 + LDVT2*M;

      IDX = 1;
      IDXC = IDX + N;
      COLTYP = IDXC + N;
      IDXP = COLTYP + N;

      // Scale.

      ORGNRM = max( ( ALPHA ).abs(), ( BETA ).abs() );
      D[NL+1] = ZERO;
      for (I = 1; I <= N; I++) { // 10
         if ( ( D( I ) ).abs() > ORGNRM ) {
            ORGNRM = ( D( I ) ).abs();
         }
      } // 10
      slascl('G', 0, 0, ORGNRM, ONE, N, 1, D, N, INFO );
      ALPHA = ALPHA / ORGNRM;
      BETA = BETA / ORGNRM;

      // Deflate singular values.

      slasd2(NL, NR, SQRE, K, D, WORK( IZ ), ALPHA, BETA, U, LDU, VT, LDVT, WORK( ISIGMA ), WORK( IU2 ), LDU2, WORK( IVT2 ), LDVT2, IWORK( IDXP ), IWORK( IDX ), IWORK( IDXC ), IDXQ, IWORK( COLTYP ), INFO );

      // Solve Secular Equation and update singular vectors.

      LDQ = K;
      slasd3(NL, NR, SQRE, K, D, WORK( IQ ), LDQ, WORK( ISIGMA ), U, LDU, WORK( IU2 ), LDU2, VT, LDVT, WORK( IVT2 ), LDVT2, IWORK( IDXC ), IWORK( COLTYP ), WORK( IZ ), INFO );

      // Report the possible convergence failure.

      if ( INFO != 0 ) {
         return;
      }

      // Unscale.

      slascl('G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO );

      // Prepare the IDXQ sorting permutation.

      N1 = K;
      N2 = N - K;
      slamrg(N1, N2, D, 1, -1, IDXQ );

      }
