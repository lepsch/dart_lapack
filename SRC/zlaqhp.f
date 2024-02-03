      SUBROUTINE ZLAQHP( UPLO, N, AP, S, SCOND, AMAX, EQUED )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, UPLO;
      int                N;
      double             AMAX, SCOND;
      // ..
      // .. Array Arguments ..
      double             S( * );
      COMPLEX*16         AP( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, THRESH;
      const              ONE = 1.0D+0, THRESH = 0.1D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, JC;
      double             CJ, LARGE, SMALL;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH;
      // EXTERNAL LSAME, DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N.LE.0 ) {
         EQUED = 'N'
         RETURN
      }

      // Initialize LARGE and SMALL.

      SMALL = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      LARGE = ONE / SMALL

      if ( SCOND.GE.THRESH .AND. AMAX.GE.SMALL .AND. AMAX.LE.LARGE ) {

         // No equilibration

         EQUED = 'N'
      } else {

         // Replace A by diag(S) * A * diag(S).

         if ( LSAME( UPLO, 'U' ) ) {

            // Upper triangle of A is stored.

            JC = 1
            for (J = 1; J <= N; J++) { // 20
               CJ = S( J )
               DO 10 I = 1, J - 1
                  AP( JC+I-1 ) = CJ*S( I )*AP( JC+I-1 )
               } // 10
               AP( JC+J-1 ) = CJ*CJ*DBLE( AP( JC+J-1 ) )
               JC = JC + J
            } // 20
         } else {

            // Lower triangle of A is stored.

            JC = 1
            for (J = 1; J <= N; J++) { // 40
               CJ = S( J )
               AP( JC ) = CJ*CJ*DBLE( AP( JC ) )
               DO 30 I = J + 1, N
                  AP( JC+I-J ) = CJ*S( I )*AP( JC+I-J )
               } // 30
               JC = JC + N - J + 1
            } // 40
         }
         EQUED = 'Y'
      }

      RETURN

      // End of ZLAQHP

      }
