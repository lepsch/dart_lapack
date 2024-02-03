      SUBROUTINE SLASD1( NL, NR, SQRE, D, ALPHA, BETA, U, LDU, VT, LDVT, IDXQ, IWORK, WORK, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDU, LDVT, NL, NR, SQRE;
      REAL               ALPHA, BETA
      // ..
      // .. Array Arguments ..
      int                IDXQ( * ), IWORK( * );
      REAL               D( * ), U( LDU, * ), VT( LDVT, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..

      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                COLTYP, I, IDX, IDXC, IDXP, IQ, ISIGMA, IU2, IVT2, IZ, K, LDQ, LDU2, LDVT2, M, N, N1, N2;
      REAL               ORGNRM
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLAMRG, SLASCL, SLASD2, SLASD3, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0

      if ( NL.LT.1 ) {
         INFO = -1
      } else if ( NR.LT.1 ) {
         INFO = -2
      } else if ( ( SQRE.LT.0 ) .OR. ( SQRE.GT.1 ) ) {
         INFO = -3
      }
      if ( INFO != 0 ) {
         xerbla('SLASD1', -INFO );
         RETURN
      }

      N = NL + NR + 1
      M = N + SQRE

      // The following values are for bookkeeping purposes only.  They are
      // integer pointers which indicate the portion of the workspace
      // used by a particular array in SLASD2 and SLASD3.

      LDU2 = N
      LDVT2 = M

      IZ = 1
      ISIGMA = IZ + M
      IU2 = ISIGMA + N
      IVT2 = IU2 + LDU2*N
      IQ = IVT2 + LDVT2*M

      IDX = 1
      IDXC = IDX + N
      COLTYP = IDXC + N
      IDXP = COLTYP + N

      // Scale.

      ORGNRM = MAX( ABS( ALPHA ), ABS( BETA ) )
      D( NL+1 ) = ZERO
      for (I = 1; I <= N; I++) { // 10
         if ( ABS( D( I ) ).GT.ORGNRM ) {
            ORGNRM = ABS( D( I ) )
         }
      } // 10
      slascl('G', 0, 0, ORGNRM, ONE, N, 1, D, N, INFO );
      ALPHA = ALPHA / ORGNRM
      BETA = BETA / ORGNRM

      // Deflate singular values.

      slasd2(NL, NR, SQRE, K, D, WORK( IZ ), ALPHA, BETA, U, LDU, VT, LDVT, WORK( ISIGMA ), WORK( IU2 ), LDU2, WORK( IVT2 ), LDVT2, IWORK( IDXP ), IWORK( IDX ), IWORK( IDXC ), IDXQ, IWORK( COLTYP ), INFO );

      // Solve Secular Equation and update singular vectors.

      LDQ = K
      slasd3(NL, NR, SQRE, K, D, WORK( IQ ), LDQ, WORK( ISIGMA ), U, LDU, WORK( IU2 ), LDU2, VT, LDVT, WORK( IVT2 ), LDVT2, IWORK( IDXC ), IWORK( COLTYP ), WORK( IZ ), INFO );

      // Report the possible convergence failure.

      if ( INFO != 0 ) {
         RETURN
      }

      // Unscale.

      slascl('G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO );

      // Prepare the IDXQ sorting permutation.

      N1 = K
      N2 = N - K
      slamrg(N1, N2, D, 1, -1, IDXQ );

      RETURN

      // End of SLASD1

      }
