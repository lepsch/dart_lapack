      void clatdf(IJOB, N, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, JPIV ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                IJOB, LDZ, N;
      double               RDSCAL, RDSUM;
      int                IPIV( * ), JPIV( * );
      Complex            RHS( * ), Z( LDZ, * );
      // ..

      int                MAXDIM;
      const              MAXDIM = 2 ;
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      int                I, INFO, J, K;
      double               RTEMP, SCALE, SMINU, SPLUS;
      Complex            BM, BP, PMONE, TEMP;
      double               RWORK( MAXDIM );
      Complex            WORK( 4*MAXDIM ), XM( MAXDIM ), XP( MAXDIM );
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CGECON, CGESC2, CLASSQ, CLASWP, CSCAL
      // ..
      // .. External Functions ..
      //- REAL               SCASUM;
      //- COMPLEX            CDOTC;
      // EXTERNAL SCASUM, CDOTC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, SQRT

      if ( IJOB != 2 ) {

         // Apply permutations IPIV to RHS

         claswp(1, RHS, LDZ, 1, N-1, IPIV, 1 );

         // Solve for L-part choosing RHS either to +1 or -1.

         PMONE = -CONE;
         for (J = 1; J <= N - 1; J++) { // 10
            BP = RHS( J ) + CONE;
            BM = RHS( J ) - CONE;
            SPLUS = ONE;

            // Look-ahead for L- part RHS(1:N-1) = +-1
            // SPLUS and SMIN computed more efficiently than in BSOLVE[1].

            SPLUS = SPLUS + double( CDOTC( N-J, Z( J+1, J ), 1, Z( J+1, J ), 1 ) );
            SMINU = double( CDOTC( N-J, Z( J+1, J ), 1, RHS( J+1 ), 1 ) );
            SPLUS = SPLUS*double( RHS( J ) );
            if ( SPLUS > SMINU ) {
               RHS[J] = BP;
            } else if ( SMINU > SPLUS ) {
               RHS[J] = BM;
            } else {

               // In this case the updating sums are equal and we can
               // choose RHS(J) +1 or -1. The first time this happens we
               // choose -1, thereafter +1. This is a simple way to get
               // good estimates of matrices like Byers well-known example
               // (see [1]). (Not done in BSOLVE.)

               RHS[J] = RHS( J ) + PMONE;
               PMONE = CONE;
            }

            // Compute the remaining r.h.s.

            TEMP = -RHS( J );
            caxpy(N-J, TEMP, Z( J+1, J ), 1, RHS( J+1 ), 1 );
         } // 10

         // Solve for U- part, lockahead for RHS(N) = +-1. This is not done
         // In BSOLVE and will hopefully give us a better estimate because
         // any ill-conditioning of the original matrix is transferred to U
         // and not to L. U(N, N) is an approximation to sigma_min(LU).

         ccopy(N-1, RHS, 1, WORK, 1 );
         WORK[N] = RHS( N ) + CONE;
         RHS[N] = RHS( N ) - CONE;
         SPLUS = ZERO;
         SMINU = ZERO;
         for (I = N; I >= 1; I--) { // 30
            TEMP = CONE / Z( I, I );
            WORK[I] = WORK( I )*TEMP;
            RHS[I] = RHS( I )*TEMP;
            for (K = I + 1; K <= N; K++) { // 20
               WORK[I] = WORK( I ) - WORK( K )*( Z( I, K )*TEMP );
               RHS[I] = RHS( I ) - RHS( K )*( Z( I, K )*TEMP );
            } // 20
            SPLUS = SPLUS + ( WORK( I ) ).abs();
            SMINU = SMINU + ( RHS( I ) ).abs();
         } // 30
         if (SPLUS > SMINU) ccopy( N, WORK, 1, RHS, 1 );

         // Apply the permutations JPIV to the computed solution (RHS)

         claswp(1, RHS, LDZ, 1, N-1, JPIV, -1 );

         // Compute the sum of squares

         classq(N, RHS, 1, RDSCAL, RDSUM );
         return;
      }

      // ENTRY IJOB = 2

      // Compute approximate nullvector XM of Z

      cgecon('I', N, Z, LDZ, ONE, RTEMP, WORK, RWORK, INFO );
      ccopy(N, WORK( N+1 ), 1, XM, 1 );

      // Compute RHS

      claswp(1, XM, LDZ, 1, N-1, IPIV, -1 );
      TEMP = CONE / sqrt( CDOTC( N, XM, 1, XM, 1 ) );
      cscal(N, TEMP, XM, 1 );
      ccopy(N, XM, 1, XP, 1 );
      caxpy(N, CONE, RHS, 1, XP, 1 );
      caxpy(N, -CONE, XM, 1, RHS, 1 );
      cgesc2(N, Z, LDZ, RHS, IPIV, JPIV, SCALE );
      cgesc2(N, Z, LDZ, XP, IPIV, JPIV, SCALE );
      if( SCASUM( N, XP, 1 ) > SCASUM( N, RHS, 1 ) ) ccopy( N, XP, 1, RHS, 1 );

      // Compute the sum of squares

      classq(N, RHS, 1, RDSCAL, RDSUM );
      return;
      }
