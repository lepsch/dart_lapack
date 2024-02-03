      void slasd2(NL, NR, SQRE, K, D, Z, ALPHA, BETA, U, LDU, VT, LDVT, DSIGMA, U2, LDU2, VT2, LDVT2, IDXP, IDX, IDXC, IDXQ, COLTYP, INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, K, LDU, LDU2, LDVT, LDVT2, NL, NR, SQRE;
      REAL               ALPHA, BETA;
      // ..
      // .. Array Arguments ..
      int                COLTYP( * ), IDX( * ), IDXC( * ), IDXP( * ), IDXQ( * );
      REAL               D( * ), DSIGMA( * ), U( LDU, * ), U2( LDU2, * ), VT( LDVT, * ), VT2( LDVT2, * ), Z( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO, EIGHT;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0, EIGHT = 8.0 ;
      // ..
      // .. Local Arrays ..
      int                CTOT( 4 ), PSM( 4 );
      // ..
      // .. Local Scalars ..
      int                CT, I, IDXI, IDXJ, IDXJP, J, JP, JPREV, K2, M, N, NLP1, NLP2;
      REAL               C, EPS, HLFTOL, S, TAU, TOL, Z1;
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLAPY2;
      // EXTERNAL SLAMCH, SLAPY2
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SLACPY, SLAMRG, SLASET, SROT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;

      if ( NL < 1 ) {
         INFO = -1;
      } else if ( NR < 1 ) {
         INFO = -2;
      } else if ( ( SQRE != 1 ) && ( SQRE != 0 ) ) {
         INFO = -3;
      }

      N = NL + NR + 1;
      M = N + SQRE;

      if ( LDU < N ) {
         INFO = -10;
      } else if ( LDVT < M ) {
         INFO = -12;
      } else if ( LDU2 < N ) {
         INFO = -15;
      } else if ( LDVT2 < M ) {
         INFO = -17;
      }
      if ( INFO != 0 ) {
         xerbla('SLASD2', -INFO );
         return;
      }

      NLP1 = NL + 1;
      NLP2 = NL + 2;

      // Generate the first part of the vector Z; and move the singular
      // values in the first part of D one position backward.

      Z1 = ALPHA*VT( NLP1, NLP1 );
      Z( 1 ) = Z1;
      DO 10 I = NL, 1, -1;
         Z( I+1 ) = ALPHA*VT( I, NLP1 );
         D( I+1 ) = D( I );
         IDXQ( I+1 ) = IDXQ( I ) + 1;
      } // 10

      // Generate the second part of the vector Z.

      for (I = NLP2; I <= M; I++) { // 20
         Z( I ) = BETA*VT( I, NLP2 );
      } // 20

      // Initialize some reference arrays.

      for (I = 2; I <= NLP1; I++) { // 30
         COLTYP( I ) = 1;
      } // 30
      for (I = NLP2; I <= N; I++) { // 40
         COLTYP( I ) = 2;
      } // 40

      // Sort the singular values into increasing order

      for (I = NLP2; I <= N; I++) { // 50
         IDXQ( I ) = IDXQ( I ) + NLP1;
      } // 50

      // DSIGMA, IDXC, IDXC, and the first column of U2
      // are used as storage space.

      for (I = 2; I <= N; I++) { // 60
         DSIGMA( I ) = D( IDXQ( I ) );
         U2( I, 1 ) = Z( IDXQ( I ) );
         IDXC( I ) = COLTYP( IDXQ( I ) );
      } // 60

      slamrg(NL, NR, DSIGMA( 2 ), 1, 1, IDX( 2 ) );

      for (I = 2; I <= N; I++) { // 70
         IDXI = 1 + IDX( I );
         D( I ) = DSIGMA( IDXI );
         Z( I ) = U2( IDXI, 1 );
         COLTYP( I ) = IDXC( IDXI );
      } // 70

      // Calculate the allowable deflation tolerance

      EPS = SLAMCH( 'Epsilon' );
      TOL = max( ABS( ALPHA ), ABS( BETA ) );
      TOL = EIGHT*EPS*max( ABS( D( N ) ), TOL );

      // There are 2 kinds of deflation -- first a value in the z-vector
      // is small, second two (or more) singular values are very close
      // together (their difference is small).

      // If the value in the z-vector is small, we simply permute the
      // array so that the corresponding singular value is moved to the
      // end.

      // If two values in the D-vector are close, we perform a two-sided
      // rotation designed to make one of the corresponding z-vector
      // entries zero, and then permute the array so that the deflated
      // singular value is moved to the end.

      // If there are multiple singular values then the problem deflates.
      // Here the number of equal singular values are found.  As each equal
      // singular value is found, an elementary reflector is computed to
      // rotate the corresponding singular subspace so that the
      // corresponding components of Z are zero in this new basis.

      K = 1;
      K2 = N + 1;
      for (J = 2; J <= N; J++) { // 80
         if ( ABS( Z( J ) ) <= TOL ) {

            // Deflate due to small z component.

            K2 = K2 - 1;
            IDXP( K2 ) = J;
            COLTYP( J ) = 4;
            if (J == N) GO TO 120;
         } else {
            JPREV = J;
            GO TO 90;
         }
      } // 80
      } // 90
      J = JPREV;
      } // 100
      J = J + 1;
      if (J > N) GO TO 110;
      if ( ABS( Z( J ) ) <= TOL ) {

         // Deflate due to small z component.

         K2 = K2 - 1;
         IDXP( K2 ) = J;
         COLTYP( J ) = 4;
      } else {

         // Check if singular values are close enough to allow deflation.

         if ( ABS( D( J )-D( JPREV ) ) <= TOL ) {

            // Deflation is possible.

            S = Z( JPREV );
            C = Z( J );

            // Find sqrt(a**2+b**2) without overflow or
            // destructive underflow.

            TAU = SLAPY2( C, S );
            C = C / TAU;
            S = -S / TAU;
            Z( J ) = TAU;
            Z( JPREV ) = ZERO;

            // Apply back the Givens rotation to the left and right
            // singular vector matrices.

            IDXJP = IDXQ( IDX( JPREV )+1 );
            IDXJ = IDXQ( IDX( J )+1 );
            if ( IDXJP <= NLP1 ) {
               IDXJP = IDXJP - 1;
            }
            if ( IDXJ <= NLP1 ) {
               IDXJ = IDXJ - 1;
            }
            srot(N, U( 1, IDXJP ), 1, U( 1, IDXJ ), 1, C, S );
            srot(M, VT( IDXJP, 1 ), LDVT, VT( IDXJ, 1 ), LDVT, C, S );
            if ( COLTYP( J ) != COLTYP( JPREV ) ) {
               COLTYP( J ) = 3;
            }
            COLTYP( JPREV ) = 4;
            K2 = K2 - 1;
            IDXP( K2 ) = JPREV;
            JPREV = J;
         } else {
            K = K + 1;
            U2( K, 1 ) = Z( JPREV );
            DSIGMA( K ) = D( JPREV );
            IDXP( K ) = JPREV;
            JPREV = J;
         }
      }
      GO TO 100;
      } // 110

      // Record the last singular value.

      K = K + 1;
      U2( K, 1 ) = Z( JPREV );
      DSIGMA( K ) = D( JPREV );
      IDXP( K ) = JPREV;

      } // 120

      // Count up the total number of the various types of columns, then
      // form a permutation which positions the four column types into
      // four groups of uniform structure (although one or more of these
      // groups may be empty).

      for (J = 1; J <= 4; J++) { // 130
         CTOT( J ) = 0;
      } // 130
      for (J = 2; J <= N; J++) { // 140
         CT = COLTYP( J );
         CTOT( CT ) = CTOT( CT ) + 1;
      } // 140

      // PSM(*) = Position in SubMatrix (of types 1 through 4)

      PSM( 1 ) = 2;
      PSM( 2 ) = 2 + CTOT( 1 );
      PSM( 3 ) = PSM( 2 ) + CTOT( 2 );
      PSM( 4 ) = PSM( 3 ) + CTOT( 3 );

      // Fill out the IDXC array so that the permutation which it induces
      // will place all type-1 columns first, all type-2 columns next,
      // then all type-3's, and finally all type-4's, starting from the
      // second column. This applies similarly to the rows of VT.

      for (J = 2; J <= N; J++) { // 150
         JP = IDXP( J );
         CT = COLTYP( JP );
         IDXC( PSM( CT ) ) = J;
         PSM( CT ) = PSM( CT ) + 1;
      } // 150

      // Sort the singular values and corresponding singular vectors into
      // DSIGMA, U2, and VT2 respectively.  The singular values/vectors
      // which were not deflated go into the first K slots of DSIGMA, U2,
      // and VT2 respectively, while those which were deflated go into the
      // last N - K slots, except that the first column/row will be treated
      // separately.

      for (J = 2; J <= N; J++) { // 160
         JP = IDXP( J );
         DSIGMA( J ) = D( JP );
         IDXJ = IDXQ( IDX( IDXP( IDXC( J ) ) )+1 );
         if ( IDXJ <= NLP1 ) {
            IDXJ = IDXJ - 1;
         }
         scopy(N, U( 1, IDXJ ), 1, U2( 1, J ), 1 );
         scopy(M, VT( IDXJ, 1 ), LDVT, VT2( J, 1 ), LDVT2 );
      } // 160

      // Determine DSIGMA(1), DSIGMA(2) and Z(1)

      DSIGMA( 1 ) = ZERO;
      HLFTOL = TOL / TWO;
      if( ABS( DSIGMA( 2 ) ) <= HLFTOL ) DSIGMA( 2 ) = HLFTOL;
      if ( M > N ) {
         Z( 1 ) = SLAPY2( Z1, Z( M ) );
         if ( Z( 1 ) <= TOL ) {
            C = ONE;
            S = ZERO;
            Z( 1 ) = TOL;
         } else {
            C = Z1 / Z( 1 );
            S = Z( M ) / Z( 1 );
         }
      } else {
         if ( ABS( Z1 ) <= TOL ) {
            Z( 1 ) = TOL;
         } else {
            Z( 1 ) = Z1;
         }
      }

      // Move the rest of the updating row to Z.

      scopy(K-1, U2( 2, 1 ), 1, Z( 2 ), 1 );

      // Determine the first column of U2, the first row of VT2 and the
      // last row of VT.

      slaset('A', N, 1, ZERO, ZERO, U2, LDU2 );
      U2( NLP1, 1 ) = ONE;
      if ( M > N ) {
         for (I = 1; I <= NLP1; I++) { // 170
            VT( M, I ) = -S*VT( NLP1, I );
            VT2( 1, I ) = C*VT( NLP1, I );
         } // 170
         for (I = NLP2; I <= M; I++) { // 180
            VT2( 1, I ) = S*VT( M, I );
            VT( M, I ) = C*VT( M, I );
         } // 180
      } else {
         scopy(M, VT( NLP1, 1 ), LDVT, VT2( 1, 1 ), LDVT2 );
      }
      if ( M > N ) {
         scopy(M, VT( M, 1 ), LDVT, VT2( M, 1 ), LDVT2 );
      }

      // The deflated singular values and their corresponding vectors go
      // into the back of D, U, and V respectively.

      if ( N > K ) {
         scopy(N-K, DSIGMA( K+1 ), 1, D( K+1 ), 1 );
         slacpy('A', N, N-K, U2( 1, K+1 ), LDU2, U( 1, K+1 ), LDU );
         slacpy('A', N-K, M, VT2( K+1, 1 ), LDVT2, VT( K+1, 1 ), LDVT );
      }

      // Copy CTOT into COLTYP for referencing in SLASD3.

      for (J = 1; J <= 4; J++) { // 190
         COLTYP( J ) = CTOT( J );
      } // 190

      return;

      // End of SLASD2

      }
