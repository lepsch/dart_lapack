      SUBROUTINE SLASD7( ICOMPQ, NL, NR, SQRE, K, D, Z, ZW, VF, VFW, VL, VLW, ALPHA, BETA, DSIGMA, IDX, IDXP, IDXQ, PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, C, S, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                GIVPTR, ICOMPQ, INFO, K, LDGCOL, LDGNUM, NL, NR, SQRE;
      REAL               ALPHA, BETA, C, S
      // ..
      // .. Array Arguments ..
      int                GIVCOL( LDGCOL, * ), IDX( * ), IDXP( * ), IDXQ( * ), PERM( * );
      REAL               D( * ), DSIGMA( * ), GIVNUM( LDGNUM, * ), VF( * ), VFW( * ), VL( * ), VLW( * ), Z( * ), ZW( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO, EIGHT
      const              ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0, EIGHT = 8.0E+0 ;
      // ..
      // .. Local Scalars ..

      int                I, IDXI, IDXJ, IDXJP, J, JP, JPREV, K2, M, N, NLP1, NLP2;
      REAL               EPS, HLFTOL, TAU, TOL, Z1
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SLAMRG, SROT, XERBLA
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLAPY2
      // EXTERNAL SLAMCH, SLAPY2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      N = NL + NR + 1
      M = N + SQRE

      if ( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) {
         INFO = -1
      } else if ( NL.LT.1 ) {
         INFO = -2
      } else if ( NR.LT.1 ) {
         INFO = -3
      } else if ( ( SQRE.LT.0 ) .OR. ( SQRE.GT.1 ) ) {
         INFO = -4
      } else if ( LDGCOL.LT.N ) {
         INFO = -22
      } else if ( LDGNUM.LT.N ) {
         INFO = -24
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'SLASD7', -INFO )
         RETURN
      }

      NLP1 = NL + 1
      NLP2 = NL + 2
      if ( ICOMPQ.EQ.1 ) {
         GIVPTR = 0
      }

      // Generate the first part of the vector Z and move the singular
      // values in the first part of D one position backward.

      Z1 = ALPHA*VL( NLP1 )
      VL( NLP1 ) = ZERO
      TAU = VF( NLP1 )
      DO 10 I = NL, 1, -1
         Z( I+1 ) = ALPHA*VL( I )
         VL( I ) = ZERO
         VF( I+1 ) = VF( I )
         D( I+1 ) = D( I )
         IDXQ( I+1 ) = IDXQ( I ) + 1
   10 CONTINUE
      VF( 1 ) = TAU

      // Generate the second part of the vector Z.

      DO 20 I = NLP2, M
         Z( I ) = BETA*VF( I )
         VF( I ) = ZERO
   20 CONTINUE

      // Sort the singular values into increasing order

      DO 30 I = NLP2, N
         IDXQ( I ) = IDXQ( I ) + NLP1
   30 CONTINUE

      // DSIGMA, IDXC, IDXC, and ZW are used as storage space.

      DO 40 I = 2, N
         DSIGMA( I ) = D( IDXQ( I ) )
         ZW( I ) = Z( IDXQ( I ) )
         VFW( I ) = VF( IDXQ( I ) )
         VLW( I ) = VL( IDXQ( I ) )
   40 CONTINUE

      CALL SLAMRG( NL, NR, DSIGMA( 2 ), 1, 1, IDX( 2 ) )

      DO 50 I = 2, N
         IDXI = 1 + IDX( I )
         D( I ) = DSIGMA( IDXI )
         Z( I ) = ZW( IDXI )
         VF( I ) = VFW( IDXI )
         VL( I ) = VLW( IDXI )
   50 CONTINUE

      // Calculate the allowable deflation tolerance

      EPS = SLAMCH( 'Epsilon' )
      TOL = MAX( ABS( ALPHA ), ABS( BETA ) )
      TOL = EIGHT*EIGHT*EPS*MAX( ABS( D( N ) ), TOL )

      // There are 2 kinds of deflation -- first a value in the z-vector
      // is small, second two (or more) singular values are very close
     t // ogether (their difference is small).

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

      K = 1
      K2 = N + 1
      DO 60 J = 2, N
         if ( ABS( Z( J ) ).LE.TOL ) {

            // Deflate due to small z component.

            K2 = K2 - 1
            IDXP( K2 ) = J
            IF( J.EQ.N ) GO TO 100
         } else {
            JPREV = J
            GO TO 70
         }
   60 CONTINUE
   70 CONTINUE
      J = JPREV
   80 CONTINUE
      J = J + 1
      IF( J.GT.N ) GO TO 90
      if ( ABS( Z( J ) ).LE.TOL ) {

         // Deflate due to small z component.

         K2 = K2 - 1
         IDXP( K2 ) = J
      } else {

         // Check if singular values are close enough to allow deflation.

         if ( ABS( D( J )-D( JPREV ) ).LE.TOL ) {

            // Deflation is possible.

            S = Z( JPREV )
            C = Z( J )

            // Find sqrt(a**2+b**2) without overflow or
            // destructive underflow.

            TAU = SLAPY2( C, S )
            Z( J ) = TAU
            Z( JPREV ) = ZERO
            C = C / TAU
            S = -S / TAU

            // Record the appropriate Givens rotation

            if ( ICOMPQ.EQ.1 ) {
               GIVPTR = GIVPTR + 1
               IDXJP = IDXQ( IDX( JPREV )+1 )
               IDXJ = IDXQ( IDX( J )+1 )
               if ( IDXJP.LE.NLP1 ) {
                  IDXJP = IDXJP - 1
               }
               if ( IDXJ.LE.NLP1 ) {
                  IDXJ = IDXJ - 1
               }
               GIVCOL( GIVPTR, 2 ) = IDXJP
               GIVCOL( GIVPTR, 1 ) = IDXJ
               GIVNUM( GIVPTR, 2 ) = C
               GIVNUM( GIVPTR, 1 ) = S
            }
            CALL SROT( 1, VF( JPREV ), 1, VF( J ), 1, C, S )
            CALL SROT( 1, VL( JPREV ), 1, VL( J ), 1, C, S )
            K2 = K2 - 1
            IDXP( K2 ) = JPREV
            JPREV = J
         } else {
            K = K + 1
            ZW( K ) = Z( JPREV )
            DSIGMA( K ) = D( JPREV )
            IDXP( K ) = JPREV
            JPREV = J
         }
      }
      GO TO 80
   90 CONTINUE

      // Record the last singular value.

      K = K + 1
      ZW( K ) = Z( JPREV )
      DSIGMA( K ) = D( JPREV )
      IDXP( K ) = JPREV

  100 CONTINUE

      // Sort the singular values into DSIGMA. The singular values which
      // were not deflated go into the first K slots of DSIGMA, except
     t // hat DSIGMA(1) is treated separately.

      DO 110 J = 2, N
         JP = IDXP( J )
         DSIGMA( J ) = D( JP )
         VFW( J ) = VF( JP )
         VLW( J ) = VL( JP )
  110 CONTINUE
      if ( ICOMPQ.EQ.1 ) {
         DO 120 J = 2, N
            JP = IDXP( J )
            PERM( J ) = IDXQ( IDX( JP )+1 )
            if ( PERM( J ).LE.NLP1 ) {
               PERM( J ) = PERM( J ) - 1
            }
  120    CONTINUE
      }

      // The deflated singular values go back into the last N - K slots of
      // D.

      CALL SCOPY( N-K, DSIGMA( K+1 ), 1, D( K+1 ), 1 )

      // Determine DSIGMA(1), DSIGMA(2), Z(1), VF(1), VL(1), VF(M), and
      // VL(M).

      DSIGMA( 1 ) = ZERO
      HLFTOL = TOL / TWO
      IF( ABS( DSIGMA( 2 ) ).LE.HLFTOL ) DSIGMA( 2 ) = HLFTOL
      if ( M.GT.N ) {
         Z( 1 ) = SLAPY2( Z1, Z( M ) )
         if ( Z( 1 ).LE.TOL ) {
            C = ONE
            S = ZERO
            Z( 1 ) = TOL
         } else {
            C = Z1 / Z( 1 )
            S = -Z( M ) / Z( 1 )
         }
         CALL SROT( 1, VF( M ), 1, VF( 1 ), 1, C, S )
         CALL SROT( 1, VL( M ), 1, VL( 1 ), 1, C, S )
      } else {
         if ( ABS( Z1 ).LE.TOL ) {
            Z( 1 ) = TOL
         } else {
            Z( 1 ) = Z1
         }
      }

      // Restore Z, VF, and VL.

      CALL SCOPY( K-1, ZW( 2 ), 1, Z( 2 ), 1 )
      CALL SCOPY( N-1, VFW( 2 ), 1, VF( 2 ), 1 )
      CALL SCOPY( N-1, VLW( 2 ), 1, VL( 2 ), 1 )

      RETURN

      // End of SLASD7

      }
