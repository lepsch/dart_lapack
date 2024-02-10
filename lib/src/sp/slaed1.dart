      void slaed1(N, D, Q, LDQ, INDXQ, RHO, CUTPNT, WORK, IWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                CUTPNT, INFO, LDQ, N;
      double               RHO;
      int                INDXQ( * ), IWORK( * );
      double               D( * ), Q( LDQ, * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                COLTYP, CPP1, I, IDLMDA, INDX, INDXC, INDXP, IQ2, IS, IW, IZ, K, N1, N2;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SLAED2, SLAED3, SLAMRG, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input parameters.

      INFO = 0;

      if ( N < 0 ) {
         INFO = -1;
      } else if ( LDQ < max( 1, N ) ) {
         INFO = -4;
      } else if ( min( 1, N / 2 ) > CUTPNT || ( N / 2 ) < CUTPNT ) {
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('SLAED1', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // The following values are integer pointers which indicate
      // the portion of the workspace
      // used by a particular array in SLAED2 and SLAED3.

      IZ = 1;
      IDLMDA = IZ + N;
      IW = IDLMDA + N;
      IQ2 = IW + N;

      INDX = 1;
      INDXC = INDX + N;
      COLTYP = INDXC + N;
      INDXP = COLTYP + N;


      // Form the z-vector which consists of the last row of Q_1 and the
      // first row of Q_2.

      scopy(CUTPNT, Q( CUTPNT, 1 ), LDQ, WORK( IZ ), 1 );
      CPP1 = CUTPNT + 1;
      scopy(N-CUTPNT, Q( CPP1, CPP1 ), LDQ, WORK( IZ+CUTPNT ), 1 );

      // Deflate eigenvalues.

      slaed2(K, N, CUTPNT, D, Q, LDQ, INDXQ, RHO, WORK( IZ ), WORK( IDLMDA ), WORK( IW ), WORK( IQ2 ), IWORK( INDX ), IWORK( INDXC ), IWORK( INDXP ), IWORK( COLTYP ), INFO );

      if (INFO != 0) GO TO 20;

      // Solve Secular Equation.

      if ( K != 0 ) {
         IS = ( IWORK( COLTYP )+IWORK( COLTYP+1 ) )*CUTPNT + ( IWORK( COLTYP+1 )+IWORK( COLTYP+2 ) )*( N-CUTPNT ) + IQ2;
         slaed3(K, N, CUTPNT, D, Q, LDQ, RHO, WORK( IDLMDA ), WORK( IQ2 ), IWORK( INDXC ), IWORK( COLTYP ), WORK( IW ), WORK( IS ), INFO );
         if (INFO != 0) GO TO 20;

      // Prepare the INDXQ sorting permutation.

         N1 = K;
         N2 = N - K;
         slamrg(N1, N2, D, 1, -1, INDXQ );
      } else {
         for (I = 1; I <= N; I++) { // 10
            INDXQ[I] = I;
         } // 10
      }

      } // 20
      }
