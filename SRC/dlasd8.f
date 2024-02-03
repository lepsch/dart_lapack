      SUBROUTINE DLASD8( ICOMPQ, K, D, Z, VF, VL, DIFL, DIFR, LDDIFR, DSIGMA, WORK, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                ICOMPQ, INFO, K, LDDIFR;
      // ..
      // .. Array Arguments ..
      double             D( * ), DIFL( * ), DIFR( LDDIFR, * ), DSIGMA( * ), VF( * ), VL( * ), WORK( * ), Z( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, IWK1, IWK2, IWK2I, IWK3, IWK3I, J;
      double             DIFLJ, DIFRJ, DJ, DSIGJ, DSIGJP, RHO, TEMP;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DLASCL, DLASD4, DLASET, XERBLA
      // ..
      // .. External Functions ..
      double             DDOT, DLAMC3, DNRM2;
      // EXTERNAL DDOT, DLAMC3, DNRM2
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
      if ( INFO.NE.0 ) {
         xerbla('DLASD8', -INFO );
         RETURN
      }

      // Quick return if possible

      if ( K.EQ.1 ) {
         D( 1 ) = ABS( Z( 1 ) )
         DIFL( 1 ) = D( 1 )
         if ( ICOMPQ.EQ.1 ) {
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

      RHO = DNRM2( K, Z, 1 )
      dlascl('G', 0, 0, RHO, ONE, K, 1, Z, K, INFO );
      RHO = RHO*RHO

      // Initialize WORK(IWK3).

      dlaset('A', K, 1, ONE, ONE, WORK( IWK3 ), K );

      // Compute the updated singular values, the arrays DIFL, DIFR,
      // and the updated Z.

      DO 40 J = 1, K
         dlasd4(K, J, DSIGMA, Z, WORK( IWK1 ), RHO, D( J ), WORK( IWK2 ), INFO );

         // If the root finder fails, report the convergence failure.

         if ( INFO.NE.0 ) {
            RETURN
         }
         WORK( IWK3I+J ) = WORK( IWK3I+J )*WORK( J )*WORK( IWK2I+J )
         DIFL( J ) = -WORK( J )
         DIFR( J, 1 ) = -WORK( J+1 )
         DO 20 I = 1, J - 1
            WORK( IWK3I+I ) = WORK( IWK3I+I )*WORK( I )* WORK( IWK2I+I ) / ( DSIGMA( I )- DSIGMA( J ) ) / ( DSIGMA( I )+ DSIGMA( J ) )
   20    CONTINUE
         DO 30 I = J + 1, K
            WORK( IWK3I+I ) = WORK( IWK3I+I )*WORK( I )* WORK( IWK2I+I ) / ( DSIGMA( I )- DSIGMA( J ) ) / ( DSIGMA( I )+ DSIGMA( J ) )
   30    CONTINUE
   40 CONTINUE

      // Compute updated Z.

      DO 50 I = 1, K
         Z( I ) = SIGN( SQRT( ABS( WORK( IWK3I+I ) ) ), Z( I ) )
   50 CONTINUE

      // Update VF and VL.

      DO 80 J = 1, K
         DIFLJ = DIFL( J )
         DJ = D( J )
         DSIGJ = -DSIGMA( J )
         if ( J.LT.K ) {
            DIFRJ = -DIFR( J, 1 )
            DSIGJP = -DSIGMA( J+1 )
         }
         WORK( J ) = -Z( J ) / DIFLJ / ( DSIGMA( J )+DJ )

         // Use calls to the subroutine DLAMC3 to enforce the parentheses
         // (x+y)+z. The goal is to prevent optimizing compilers
         // from doing x+(y+z).

         DO 60 I = 1, J - 1
            WORK( I ) = Z( I ) / ( DLAMC3( DSIGMA( I ), DSIGJ )-DIFLJ ) / ( DSIGMA( I )+DJ )
   60    CONTINUE
         DO 70 I = J + 1, K
            WORK( I ) = Z( I ) / ( DLAMC3( DSIGMA( I ), DSIGJP )+DIFRJ ) / ( DSIGMA( I )+DJ )
   70    CONTINUE
         TEMP = DNRM2( K, WORK, 1 )
         WORK( IWK2I+J ) = DDOT( K, WORK, 1, VF, 1 ) / TEMP
         WORK( IWK3I+J ) = DDOT( K, WORK, 1, VL, 1 ) / TEMP
         if ( ICOMPQ.EQ.1 ) {
            DIFR( J, 2 ) = TEMP
         }
   80 CONTINUE

      dcopy(K, WORK( IWK2 ), 1, VF, 1 );
      dcopy(K, WORK( IWK3 ), 1, VL, 1 );

      RETURN

      // End of DLASD8

      }
