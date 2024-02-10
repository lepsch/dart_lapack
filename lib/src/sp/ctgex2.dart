      void ctgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, J1, Box<int> INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               WANTQ, WANTZ;
      int                INFO, J1, LDA, LDB, LDQ, LDZ, N;
      Complex            A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * );
      // ..

      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      double               TWENTY;
      const              TWENTY = 2.0e+1 ;
      int                LDST;
      const              LDST = 2 ;
      bool               WANDS;
      const              WANDS = true ;
      bool               STRONG, WEAK;
      int                I, M;
      double               CQ, CZ, EPS, SA, SB, SCALE, SMLNUM, SUM, THRESHA, THRESHB;
      Complex            CDUM, F, G, SQ, SZ;
      Complex            S( LDST, LDST ), T( LDST, LDST ), WORK( 8 );
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACPY, CLARTG, CLASSQ, CROT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, MAX, REAL, SQRT

      INFO = 0;

      // Quick return if possible

      if (N <= 1) return;

      M = LDST;
      WEAK = false;
      STRONG = false;

      // Make a local copy of selected block in (A, B)

      clacpy('Full', M, M, A( J1, J1 ), LDA, S, LDST );
      clacpy('Full', M, M, B( J1, J1 ), LDB, T, LDST );

      // Compute the threshold for testing the acceptance of swapping.

      EPS = SLAMCH( 'P' );
      SMLNUM = SLAMCH( 'S' ) / EPS;
      SCALE = double( CZERO );
      SUM = double( CONE );
      clacpy('Full', M, M, S, LDST, WORK, M );
      clacpy('Full', M, M, T, LDST, WORK( M*M+1 ), M );
      classq(M*M, WORK, 1, SCALE, SUM );
      SA = SCALE*sqrt( SUM );
      SCALE = CZERO.toDouble();
      SUM = CONE.toDouble();
      classq(M*M, WORK(M*M+1), 1, SCALE, SUM );
      SB = SCALE*sqrt( SUM );

      // THRES has been changed from
      //    THRESH = max( TEN*EPS*SA, SMLNUM )
      // to
      //    THRESH = max( TWENTY*EPS*SA, SMLNUM )
      // on 04/01/10.
      // "Bug" reported by Ondra Kamenik, confirmed by Julie Langou, fixed by
      // Jim Demmel and Guillaume Revy. See forum post 1783.

      THRESHA = max( TWENTY*EPS*SA, SMLNUM );
      THRESHB = max( TWENTY*EPS*SB, SMLNUM );

      // Compute unitary QL and RQ that swap 1-by-1 and 1-by-1 blocks
      // using Givens rotations and perform the swap tentatively.

      F = S( 2, 2 )*T( 1, 1 ) - T( 2, 2 )*S( 1, 1 );
      G = S( 2, 2 )*T( 1, 2 ) - T( 2, 2 )*S( 1, 2 );
      SA = ( S( 2, 2 ) ).abs() * ( T( 1, 1 ) ).abs();
      SB = ( S( 1, 1 ) ).abs() * ( T( 2, 2 ) ).abs();
      clartg(G, F, CZ, SZ, CDUM );
      SZ = -SZ;
      crot(2, S( 1, 1 ), 1, S( 1, 2 ), 1, CZ, CONJG( SZ ) );
      crot(2, T( 1, 1 ), 1, T( 1, 2 ), 1, CZ, CONJG( SZ ) );
      if ( SA >= SB ) {
         clartg(S( 1, 1 ), S( 2, 1 ), CQ, SQ, CDUM );
      } else {
         clartg(T( 1, 1 ), T( 2, 1 ), CQ, SQ, CDUM );
      }
      crot(2, S( 1, 1 ), LDST, S( 2, 1 ), LDST, CQ, SQ );
      crot(2, T( 1, 1 ), LDST, T( 2, 1 ), LDST, CQ, SQ );

      // Weak stability test: |S21| <= O(EPS F-norm((A)))
      //                      and  |T21| <= O(EPS F-norm((B)))

      WEAK = ( S( 2, 1 ) ).abs() <= THRESHA && ( T( 2, 1 ) ).abs() <= THRESHB       IF( !WEAK ) GO TO 20;

      if ( WANDS ) {

         // Strong stability test:
         //    F-norm((A-QL**H*S*QR, B-QL**H*T*QR)) <= O(EPS*F-norm((A, B)))

         clacpy('Full', M, M, S, LDST, WORK, M );
         clacpy('Full', M, M, T, LDST, WORK( M*M+1 ), M );
         crot(2, WORK, 1, WORK( 3 ), 1, CZ, -CONJG( SZ ) );
         crot(2, WORK( 5 ), 1, WORK( 7 ), 1, CZ, -CONJG( SZ ) );
         crot(2, WORK, 2, WORK( 2 ), 2, CQ, -SQ );
         crot(2, WORK( 5 ), 2, WORK( 6 ), 2, CQ, -SQ );
         for (I = 1; I <= 2; I++) { // 10
            WORK[I] = WORK( I ) - A( J1+I-1, J1 );
            WORK[I+2] = WORK( I+2 ) - A( J1+I-1, J1+1 );
            WORK[I+4] = WORK( I+4 ) - B( J1+I-1, J1 );
            WORK[I+6] = WORK( I+6 ) - B( J1+I-1, J1+1 );
         } // 10
         SCALE = CZERO.toDouble();
         SUM = CONE.toDouble();
         classq(M*M, WORK, 1, SCALE, SUM );
         SA = SCALE*sqrt( SUM );
         SCALE = CZERO.toDouble();
         SUM = CONE.toDouble();
         classq(M*M, WORK(M*M+1), 1, SCALE, SUM );
         SB = SCALE*sqrt( SUM );
         STRONG = SA <= THRESHA && SB <= THRESHB;
         if ( !STRONG) GO TO 20;
      }

      // If the swap is accepted ("weakly" and "strongly"), apply the
      // equivalence transformations to the original matrix pair (A,B)

      crot(J1+1, A( 1, J1 ), 1, A( 1, J1+1 ), 1, CZ, CONJG( SZ ) );
      crot(J1+1, B( 1, J1 ), 1, B( 1, J1+1 ), 1, CZ, CONJG( SZ ) );
      crot(N-J1+1, A( J1, J1 ), LDA, A( J1+1, J1 ), LDA, CQ, SQ );
      crot(N-J1+1, B( J1, J1 ), LDB, B( J1+1, J1 ), LDB, CQ, SQ );

      // Set  N1 by N2 (2,1) blocks to 0

      A[J1+1][J1] = CZERO;
      B[J1+1][J1] = CZERO;

      // Accumulate transformations into Q and Z if requested.

      if (WANTZ) crot( N, Z( 1, J1 ), 1, Z( 1, J1+1 ), 1, CZ, CONJG( SZ ) );
      IF( WANTQ ) crot( N, Q( 1, J1 ), 1, Q( 1, J1+1 ), 1, CQ, CONJG( SQ ) );

      // Exit with INFO = 0 if swap was successfully performed.

      return;

      // Exit with INFO = 1 if swap was rejected.

      } // 20
      INFO = 1;
      }
