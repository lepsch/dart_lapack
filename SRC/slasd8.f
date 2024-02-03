      SUBROUTINE SLASD8( ICOMPQ, K, D, Z, VF, VL, DIFL, DIFR, LDDIFR, DSIGMA, WORK, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                ICOMPQ, INFO, K, LDDIFR;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), DIFL( * ), DIFR( LDDIFR, * ), DSIGMA( * ), VF( * ), VL( * ), WORK( * ), Z( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, IWK1, IWK2, IWK2I, IWK3, IWK3I, J;
      REAL               DIFLJ, DIFRJ, DJ, DSIGJ, DSIGJP, RHO, TEMP
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SLASCL, SLASD4, SLASET, XERBLA
      // ..
      // .. External Functions ..
      REAL               SDOT, SLAMC3, SNRM2
      // EXTERNAL SDOT, SLAMC3, SNRM2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SIGN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0

      if ( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) {
         INFO = -1
      } else if ( K.LT.1 ) {
         INFO = -2
      } else if ( LDDIFR.LT.K ) {
         INFO = -9
      }
      if ( INFO != 0 ) {
         xerbla('SLASD8', -INFO );
         RETURN
      }

      // Quick return if possible

      if ( K == 1 ) {
         D( 1 ) = ABS( Z( 1 ) )
         DIFL( 1 ) = D( 1 )
         if ( ICOMPQ == 1 ) {
            DIFL( 2 ) = ONE
            DIFR( 1, 2 ) = ONE
         }
         RETURN
      }

      // Book keeping.

      IWK1 = 1
      IWK2 = IWK1 + K
      IWK3 = IWK2 + K
      IWK2I = IWK2 - 1
      IWK3I = IWK3 - 1

      // Normalize Z.

      RHO = SNRM2( K, Z, 1 )
      slascl('G', 0, 0, RHO, ONE, K, 1, Z, K, INFO );
      RHO = RHO*RHO

      // Initialize WORK(IWK3).

      slaset('A', K, 1, ONE, ONE, WORK( IWK3 ), K );

      // Compute the updated singular values, the arrays DIFL, DIFR,
      // and the updated Z.

      for (J = 1; J <= K; J++) { // 40
         slasd4(K, J, DSIGMA, Z, WORK( IWK1 ), RHO, D( J ), WORK( IWK2 ), INFO );

         // If the root finder fails, report the convergence failure.

         if ( INFO != 0 ) {
            RETURN
         }
         WORK( IWK3I+J ) = WORK( IWK3I+J )*WORK( J )*WORK( IWK2I+J )
         DIFL( J ) = -WORK( J )
         DIFR( J, 1 ) = -WORK( J+1 )
         for (I = 1; I <= J - 1; I++) { // 20
            WORK( IWK3I+I ) = WORK( IWK3I+I )*WORK( I )* WORK( IWK2I+I ) / ( DSIGMA( I )- DSIGMA( J ) ) / ( DSIGMA( I )+ DSIGMA( J ) )
         } // 20
         for (I = J + 1; I <= K; I++) { // 30
            WORK( IWK3I+I ) = WORK( IWK3I+I )*WORK( I )* WORK( IWK2I+I ) / ( DSIGMA( I )- DSIGMA( J ) ) / ( DSIGMA( I )+ DSIGMA( J ) )
         } // 30
      } // 40

      // Compute updated Z.

      for (I = 1; I <= K; I++) { // 50
         Z( I ) = SIGN( SQRT( ABS( WORK( IWK3I+I ) ) ), Z( I ) )
      } // 50

      // Update VF and VL.

      for (J = 1; J <= K; J++) { // 80
         DIFLJ = DIFL( J )
         DJ = D( J )
         DSIGJ = -DSIGMA( J )
         if ( J.LT.K ) {
            DIFRJ = -DIFR( J, 1 )
            DSIGJP = -DSIGMA( J+1 )
         }
         WORK( J ) = -Z( J ) / DIFLJ / ( DSIGMA( J )+DJ )

         // Use calls to the subroutine SLAMC3 to enforce the parentheses
         // (x+y)+z. The goal is to prevent optimizing compilers
         // from doing x+(y+z).

         for (I = 1; I <= J - 1; I++) { // 60
            WORK( I ) = Z( I ) / ( SLAMC3( DSIGMA( I ), DSIGJ )-DIFLJ ) / ( DSIGMA( I )+DJ )
         } // 60
         for (I = J + 1; I <= K; I++) { // 70
            WORK( I ) = Z( I ) / ( SLAMC3( DSIGMA( I ), DSIGJP )+DIFRJ ) / ( DSIGMA( I )+DJ )
         } // 70
         TEMP = SNRM2( K, WORK, 1 )
         WORK( IWK2I+J ) = SDOT( K, WORK, 1, VF, 1 ) / TEMP
         WORK( IWK3I+J ) = SDOT( K, WORK, 1, VL, 1 ) / TEMP
         if ( ICOMPQ == 1 ) {
            DIFR( J, 2 ) = TEMP
         }
      } // 80

      scopy(K, WORK( IWK2 ), 1, VF, 1 );
      scopy(K, WORK( IWK3 ), 1, VL, 1 );

      RETURN

      // End of SLASD8

      }
