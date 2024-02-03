      void slasd6(ICOMPQ, NL, NR, SQRE, D, VF, VL, ALPHA, BETA, IDXQ, PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, POLES, DIFL, DIFR, Z, K, C, S, WORK, IWORK, INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                GIVPTR, ICOMPQ, INFO, K, LDGCOL, LDGNUM, NL, NR, SQRE;
      REAL               ALPHA, BETA, C, S;
      // ..
      // .. Array Arguments ..
      int                GIVCOL( LDGCOL, * ), IDXQ( * ), IWORK( * ), PERM( * );
      REAL               D( * ), DIFL( * ), DIFR( * ), GIVNUM( LDGNUM, * ), POLES( LDGNUM, * ), VF( * ), VL( * ), WORK( * ), Z( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I, IDX, IDXC, IDXP, ISIGMA, IVFW, IVLW, IW, M, N, N1, N2;
      REAL               ORGNRM;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SLAMRG, SLASCL, SLASD7, SLASD8, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      N = NL + NR + 1;
      M = N + SQRE;

      if ( ( ICOMPQ < 0 ) || ( ICOMPQ > 1 ) ) {
         INFO = -1;
      } else if ( NL < 1 ) {
         INFO = -2;
      } else if ( NR < 1 ) {
         INFO = -3;
      } else if ( ( SQRE < 0 ) || ( SQRE > 1 ) ) {
         INFO = -4;
      } else if ( LDGCOL < N ) {
         INFO = -14;
      } else if ( LDGNUM < N ) {
         INFO = -16;
      }
      if ( INFO != 0 ) {
         xerbla('SLASD6', -INFO );
         return;
      }

      // The following values are for bookkeeping purposes only.  They are
      // integer pointers which indicate the portion of the workspace
      // used by a particular array in SLASD7 and SLASD8.

      ISIGMA = 1;
      IW = ISIGMA + N;
      IVFW = IW + M;
      IVLW = IVFW + M;

      IDX = 1;
      IDXC = IDX + N;
      IDXP = IDXC + N;

      // Scale.

      ORGNRM = max( ABS( ALPHA ), ABS( BETA ) );
      D( NL+1 ) = ZERO;
      for (I = 1; I <= N; I++) { // 10
         if ( ABS( D( I ) ) > ORGNRM ) {
            ORGNRM = ABS( D( I ) );
         }
      } // 10
      slascl('G', 0, 0, ORGNRM, ONE, N, 1, D, N, INFO );
      ALPHA = ALPHA / ORGNRM;
      BETA = BETA / ORGNRM;

      // Sort and Deflate singular values.

      slasd7(ICOMPQ, NL, NR, SQRE, K, D, Z, WORK( IW ), VF, WORK( IVFW ), VL, WORK( IVLW ), ALPHA, BETA, WORK( ISIGMA ), IWORK( IDX ), IWORK( IDXP ), IDXQ, PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, C, S, INFO );

      // Solve Secular Equation, compute DIFL, DIFR, and update VF, VL.

      slasd8(ICOMPQ, K, D, Z, VF, VL, DIFL, DIFR, LDGNUM, WORK( ISIGMA ), WORK( IW ), INFO );

      // Report the possible convergence failure.

      if ( INFO != 0 ) {
         return;
      }

      // Save the poles if ICOMPQ = 1.

      if ( ICOMPQ == 1 ) {
         scopy(K, D, 1, POLES( 1, 1 ), 1 );
         scopy(K, WORK( ISIGMA ), 1, POLES( 1, 2 ), 1 );
      }

      // Unscale.

      slascl('G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO );

      // Prepare the IDXQ sorting permutation.

      N1 = K;
      N2 = N - K;
      slamrg(N1, N2, D, 1, -1, IDXQ );

      return;

      // End of SLASD6

      }
