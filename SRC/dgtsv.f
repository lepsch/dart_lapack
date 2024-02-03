      SUBROUTINE DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      double             B( LDB, * ), D( * ), DL( * ), DU( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      double             FACT, TEMP;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Executable Statements ..

      INFO = 0
      if ( N.LT.0 ) {
         INFO = -1
      } else if ( NRHS.LT.0 ) {
         INFO = -2
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -7
      }
      if ( INFO.NE.0 ) {
         xerbla('DGTSV ', -INFO );
         RETURN
      }

      IF( N.EQ.0 ) RETURN

      if ( NRHS.EQ.1 ) {
         DO 10 I = 1, N - 2
            if ( ABS( D( I ) ).GE.ABS( DL( I ) ) ) {

               // No row interchange required

               if ( D( I ).NE.ZERO ) {
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  B( I+1, 1 ) = B( I+1, 1 ) - FACT*B( I, 1 )
               } else {
                  INFO = I
                  RETURN
               }
               DL( I ) = ZERO
            } else {

               // Interchange rows I and I+1

               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DL( I ) = DU( I+1 )
               DU( I+1 ) = -FACT*DL( I )
               DU( I ) = TEMP
               TEMP = B( I, 1 )
               B( I, 1 ) = B( I+1, 1 )
               B( I+1, 1 ) = TEMP - FACT*B( I+1, 1 )
            }
   10    CONTINUE
         if ( N.GT.1 ) {
            I = N - 1
            if ( ABS( D( I ) ).GE.ABS( DL( I ) ) ) {
               if ( D( I ).NE.ZERO ) {
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  B( I+1, 1 ) = B( I+1, 1 ) - FACT*B( I, 1 )
               } else {
                  INFO = I
                  RETURN
               }
            } else {
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DU( I ) = TEMP
               TEMP = B( I, 1 )
               B( I, 1 ) = B( I+1, 1 )
               B( I+1, 1 ) = TEMP - FACT*B( I+1, 1 )
            }
         }
         if ( D( N ).EQ.ZERO ) {
            INFO = N
            RETURN
         }
      } else {
         DO 40 I = 1, N - 2
            if ( ABS( D( I ) ).GE.ABS( DL( I ) ) ) {

               // No row interchange required

               if ( D( I ).NE.ZERO ) {
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  DO 20 J = 1, NRHS
                     B( I+1, J ) = B( I+1, J ) - FACT*B( I, J )
   20             CONTINUE
               } else {
                  INFO = I
                  RETURN
               }
               DL( I ) = ZERO
            } else {

               // Interchange rows I and I+1

               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DL( I ) = DU( I+1 )
               DU( I+1 ) = -FACT*DL( I )
               DU( I ) = TEMP
               DO 30 J = 1, NRHS
                  TEMP = B( I, J )
                  B( I, J ) = B( I+1, J )
                  B( I+1, J ) = TEMP - FACT*B( I+1, J )
   30          CONTINUE
            }
   40    CONTINUE
         if ( N.GT.1 ) {
            I = N - 1
            if ( ABS( D( I ) ).GE.ABS( DL( I ) ) ) {
               if ( D( I ).NE.ZERO ) {
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  DO 50 J = 1, NRHS
                     B( I+1, J ) = B( I+1, J ) - FACT*B( I, J )
   50             CONTINUE
               } else {
                  INFO = I
                  RETURN
               }
            } else {
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DU( I ) = TEMP
               DO 60 J = 1, NRHS
                  TEMP = B( I, J )
                  B( I, J ) = B( I+1, J )
                  B( I+1, J ) = TEMP - FACT*B( I+1, J )
   60          CONTINUE
            }
         }
         if ( D( N ).EQ.ZERO ) {
            INFO = N
            RETURN
         }
      }

      // Back solve with the matrix U from the factorization.

      if ( NRHS.LE.2 ) {
         J = 1
   70    CONTINUE
         B( N, J ) = B( N, J ) / D( N )
         IF( N.GT.1 ) B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / D( N-1 )
         DO 80 I = N - 2, 1, -1
            B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DL( I )* B( I+2, J ) ) / D( I )
   80    CONTINUE
         if ( J.LT.NRHS ) {
            J = J + 1
            GO TO 70
         }
      } else {
         DO 100 J = 1, NRHS
            B( N, J ) = B( N, J ) / D( N )
            IF( N.GT.1 ) B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / D( N-1 )
            DO 90 I = N - 2, 1, -1
               B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DL( I )* B( I+2, J ) ) / D( I )
   90       CONTINUE
  100    CONTINUE
      }

      RETURN

      // End of DGTSV

      }
