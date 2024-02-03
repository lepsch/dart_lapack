      SUBROUTINE ZBDT03( UPLO, N, KD, D, E, U, LDU, S, VT, LDVT, WORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                KD, LDU, LDVT, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * ), S( * );
      COMPLEX*16         U( LDU, * ), VT( LDVT, * ), WORK( * )
      // ..

* ======================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      double             BNORM, EPS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX;
      double             DLAMCH, DZASUM;
      // EXTERNAL LSAME, IDAMAX, DLAMCH, DZASUM
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      RESID = ZERO
      IF( N.LE.0 ) RETURN

      // Compute B - U * S * V' one column at a time.

      BNORM = ZERO
      if ( KD.GE.1 ) {

         // B is bidiagonal.

         if ( LSAME( UPLO, 'U' ) ) {

            // B is upper bidiagonal.

            DO 20 J = 1, N
               DO 10 I = 1, N
                  WORK( N+I ) = S( I )*VT( I, J )
   10          CONTINUE
               zgemv('No transpose', N, N, -DCMPLX( ONE ), U, LDU, WORK( N+1 ), 1, DCMPLX( ZERO ), WORK, 1 );
               WORK( J ) = WORK( J ) + D( J )
               if ( J.GT.1 ) {
                  WORK( J-1 ) = WORK( J-1 ) + E( J-1 )
                  BNORM = MAX( BNORM, ABS( D( J ) )+ABS( E( J-1 ) ) )
               } else {
                  BNORM = MAX( BNORM, ABS( D( J ) ) )
               }
               RESID = MAX( RESID, DZASUM( N, WORK, 1 ) )
   20       CONTINUE
         } else {

            // B is lower bidiagonal.

            DO 40 J = 1, N
               DO 30 I = 1, N
                  WORK( N+I ) = S( I )*VT( I, J )
   30          CONTINUE
               zgemv('No transpose', N, N, -DCMPLX( ONE ), U, LDU, WORK( N+1 ), 1, DCMPLX( ZERO ), WORK, 1 );
               WORK( J ) = WORK( J ) + D( J )
               if ( J.LT.N ) {
                  WORK( J+1 ) = WORK( J+1 ) + E( J )
                  BNORM = MAX( BNORM, ABS( D( J ) )+ABS( E( J ) ) )
               } else {
                  BNORM = MAX( BNORM, ABS( D( J ) ) )
               }
               RESID = MAX( RESID, DZASUM( N, WORK, 1 ) )
   40       CONTINUE
         }
      } else {

         // B is diagonal.

         DO 60 J = 1, N
            DO 50 I = 1, N
               WORK( N+I ) = S( I )*VT( I, J )
   50       CONTINUE
            zgemv('No transpose', N, N, -DCMPLX( ONE ), U, LDU, WORK( N+1 ), 1, DCMPLX( ZERO ), WORK, 1 );
            WORK( J ) = WORK( J ) + D( J )
            RESID = MAX( RESID, DZASUM( N, WORK, 1 ) )
   60    CONTINUE
         J = IDAMAX( N, D, 1 )
         BNORM = ABS( D( J ) )
      }

      // Compute norm(B - U * S * V') / ( n * norm(B) * EPS )

      EPS = DLAMCH( 'Precision' )

      if ( BNORM.LE.ZERO ) {
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      } else {
         if ( BNORM.GE.RESID ) {
            RESID = ( RESID / BNORM ) / ( DBLE( N )*EPS )
         } else {
            if ( BNORM.LT.ONE ) {
               RESID = ( MIN( RESID, DBLE( N )*BNORM ) / BNORM ) / ( DBLE( N )*EPS )
            } else {
               RESID = MIN( RESID / BNORM, DBLE( N ) ) / ( DBLE( N )*EPS )
            }
         }
      }

      RETURN

      // End of ZBDT03

      }
